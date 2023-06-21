library(exact2x2)

oldCrossCheck = function() {
  Tab1 = read_csv("genome_indices_from Tim_22Apr2021.csv", guess_max = 2^20)
  Tab2 = read_csv("UNITAID_FIND_WHO_SimplifiedTable_FINAL_31Mar2021.csv", guess_max = 2^20)
  Tab3 = read_csv("Genome_Indices_FINAL_22Apr2021.csv", guess_max = 2^20)
  coln = colnames(Tab2)
  firstRow = Tab2 %>% 
    slice(1)
  coln[str_sub(coln, 1, 1) == "X"] = firstRow[str_sub(coln, 1, 1) == "X"]
  colnames(Tab2) = unlist(coln)
  Tab2 = Tab2 %>%
    dplyr::filter(c(FALSE, rep(TRUE, nrow(Tab2) - 1)))
  total2 = Tab2 %>%
    select(starts_with("MUT")) %>%
    mutate_all(as.integer) %>%
    rowSums(na.rm = TRUE)
  Tab2 = Tab2 %>%
    dplyr::filter(total2 > 0)
  Tab2 = Tab2 %>%
    mutate(variant = `variant (common_name)` %>% str_remove_all("\\ .*"))
  Tab2E = Tab2 %>%
    inner_join(Tab1, by = c("variant" = "variant", "drug" = "drug"))
  Tab2E = Tab2E %>%
    select(c(drug:`Genome position`, `FINAL CONFIDENCE GRADING`:found_in_S)) %>%
    select(-c(tier, `variant (common_name)`)) %>%
    select(gene, variant, `Genome position`, ref, alt, found_in_R, found_in_S, `FINAL CONFIDENCE GRADING`, everything())
  Tab2E = Tab2E %>%
    rename(gene_name = gene, genome_index = `Genome position`, ref_nt = ref, alt_nt = alt, Total_R = found_in_R, 
           Total_S = found_in_S, Conf_Grade = `FINAL CONFIDENCE GRADING`)
  Tab2E = Tab2E %>%
    mutate(mutation = str_sub(variant, nchar(gene_name) + 2)) %>%
    mutate(codon_number = str_extract(mutation, "^[A-Z]{1}[0-9]+[A-Z\\!]{1}$") %>% str_sub(2, -2) %>% as.integer()) %>%
    mutate(ref_aa = str_extract(mutation, "^[A-Z]"), alt_aa = str_extract(mutation, "[A-Z\\!]$")) %>%
    select(-mutation)
  Tab2E = Tab2E %>%
    select(gene_name, variant, codon_number, genome_index, ref_nt, alt_nt, ref_aa, alt_aa, Total_R, Total_S, everything())
  Tab2E = Tab2E %>%
    dplyr::filter(Conf_Grade != "combo")
  Tab2F = Tab2E %>%
    group_by(variant) %>%
    group_split() %>%
    lapply(function(x) {pivot_wider(x, names_from = drug, values_from = c(Total_R, Total_S, Conf_Grade))})
  Tab2G = do.call(bind_rows, Tab2F)
  TotalR = Tab2G %>%
    select(starts_with("Total_R")) %>%
    rowSums(na.rm = TRUE)
  TotalS = Tab2G %>%
    select(starts_with("Total_S")) %>%
    rowSums(na.rm = TRUE)
  Tab2G = Tab2G %>%
    mutate(Total_R = TotalR, Total_S = TotalS) %>%
    select(gene_name:alt_aa, Total_R, Total_S, Total_R_AMI:Conf_Grade_RIF, genome_index1:alt3)
  Tab2G = Tab2G %>%
    mutate(alt_codon_number = variant %>% str_extract("\\_[0-9]+\\_(ins|del)\\_") %>% str_remove_all("(ins|del)") %>% 
             str_remove_all("\\_") %>% as.integer() %>% add(2) %>% divide_by_int(3) %>% str_c("(", ., ")"))
  Tab3E = Tab3 %>%
    dplyr::filter(is.na(ref_aa) | is.na(alt_aa) | ref_aa != alt_aa)
  Vars3 = unique(Tab3E$variant)
  Vars2 = unique(Tab2G$variant)
  missingVars = setdiff(Vars2, Vars3)
  print(paste("There are", length(missingVars), "variants that only appear in the initial catalogue"))
  print(paste("They are:",  paste(missingVars, collapse = ", ")))
  omittedVars = setdiff(Vars3, Vars2)
  print(paste("There are", length(omittedVars), "variants that only appear in the final catalogue"))
  print(paste("They are:",  paste(omittedVars, collapse = ", ")))
  Tab2H = Tab2G %>%
    dplyr::filter(variant %in% Vars3) %>%
    arrange(variant, ref_nt, alt_nt)
  ### Source: https://mycobrowser.epfl.ch/releases (release 4, tab-separated files)
  extraTab1 = read_tsv("Mycobacterium_tuberculosis_H37Rv_txt_v4.txt", guess_max = 2^20) %>%
    select(c(Locus, Name)) %>%
    set_colnames(c("gene_locus", "gene_name"))
  ### Source: https://www.ddbj.nig.ac.jp/ddbj/code-e.html
  extraTab2 = read_tsv("AminoAcids.txt") %>%
    set_colnames(c("Abbreviation", "abbr", "name")) %>%
    select(c(Abbreviation, abbr)) %>%
    slice(-27) %>%
    bind_rows(tibble(Abbreviation = "Ter", abbr = "!"))
  Tab2H = Tab2H %>%
    left_join(extraTab1, by = "gene_name") %>%
    mutate_at(c("ref_aa", "alt_aa"), ~{MM = match(., extraTab2$abbr); out = extraTab2[MM, "Abbreviation"]; out})
  Tab3E = Tab3E %>%
    arrange(variant, ref_nt, alt_nt)
  coln = colnames(Tab2H)
  coln[11:55] = paste0(str_sub(coln[11:55], -3, -1), "_", str_sub(coln[11:55], 1,-5))
  coln[11:55] = str_remove(coln[11:55], "_Total")
  coln[56:64] = paste0("detail_", coln[56:64])
  colnames(Tab2H) = coln
  for (ind in 1:nrow(Tab2H)) {
    if (ind %% 1000 == 0) { print(ind) }
    curRow = Tab2H %>%
      slice(ind)
    curInds = curRow %>%
      select(c(detail_genome_index1, detail_genome_index2, detail_genome_index3))
    curInds = curInds[!is.na(curInds)]
    if (length(curInds) > 0) {
      Tab2H$genome_index[ind] = paste(curInds, collapse = "") %>% 
        as.numeric()
    }
    curCodons = curRow %>%
      select(c(codon_number, alt_codon_number))
    curCodons = curCodons[!is.na(curCodons)]
    if (length(curCodons) > 0) {
      stopifnot(length(curCodons) == 1)
      Tab2H$codon_number[ind] = paste(curCodons, collapse = "")
    }
  }
  Tab2H = Tab2H %>%
    select(-alt_codon_number)
  coln = colnames(Tab2H)
  stopifnot(colnames(Tab3E)[61] == "detail_ref2")
  Tab3F = Tab3E %>%
    dplyr::filter(!(variant %in% Tab2H$variant)) %>%
    group_by(variant) %>%
    mutate(N = n()) %>%
    ungroup %>%
    arrange(-N)
  Tab3H = Tab3E %>%
    dplyr::filter(variant %in% Tab2H$variant)
  ### Special exception for Rv1258c_c-23t
  Tab3H = Tab3H %>%
    dplyr::filter(variant != "Rv1258c_c-23t" | Total_S != 0)
  for (Col in coln) {
    print(Col)
    curColTrain = Tab3H %>%
      select(one_of(Col))
    curColTest = Tab2H %>%
      select(one_of(Col))
    if (Col == "genome_index") {
      stopifnot(all.equal.numeric(curColTrain, curColTest))
    } else {
      n1 = sum(curColTrain != curColTest, na.rm = TRUE)
    }
    if (n1 > 0) { stop(paste("There are", n1, "discrepancies")) }
    n2 = sum(is.na(curColTrain) & !is.na(curColTest))
    if (n2 > 0) { stop(paste("There are", n2, "missing values in training set only")) }
    n3 = sum(!is.na(curColTrain) & is.na(curColTest))
    if (n3 > 0) { print(paste("There are", n3, "missing values in testing set only")) }
  }
}

processTabs = function(inputFiles = c("Full_dataset_ALL.csv", "Full_dataset_WHO.csv")) {
  outputFiles = vector("list", length(inputFiles))
  for (ind in 1:length(inputFiles)) {
    inputFile = inputFiles[ind]
    print(paste("Processing input file", ind, ":", inputFile))
    Tab = read_csv(inputFile) %>%
      select(c("drug", "tier", "variant", starts_with("...")))
    coln = colnames(Tab)
    colnT = c(coln[1:3], Tab %>% 
                slice(1) %>% 
                select(-(1:3))) 
    colnT %<>% 
      unlist %>% 
      magrittr::set_names(NULL) %>%
      str_remove_all("MUT") %>% 
      str_remove_all("present as") %>% 
      str_remove_all("pheno") %>% 
      str_remove_all(" ")
    TabT = Tab %>%
      set_colnames(colnT) %>%
      slice(-1) %>%
      rowid_to_column("origIndex")
    mutInfo = TabT %>%
      select(c(origIndex, drug, variant))
    TabT %<>%
      select(-c(drug, variant)) %>%
      mutate_all(as.integer) %>%
      inner_join(mutInfo, by = "origIndex") %>%
      mutate(SOLO_S = SOLO_SorR - SOLO_R)
    TabT %<>%
      computeCatalogueStats()
    outputFile = str_replace(inputFile, ".csv", paste0("_Processed.csv"))
    write_csv(TabT, outputFile)
    outputFiles[[ind]] = TabT
  }
  outputFiles
}

## Compute the set of mutations that are in set A but not in set B
postprocessNeutralLists = function() {
  masterTab = read_csv("allVariantsWithSetMarks.csv")
  microTab = masterTab %>% 
    dplyr::filter(setA & !setB) %>%
    group_by(variant, drug) %>% 
    slice(1) %>%
    select(variant, drug)
  ## Write that set of mutations into a separate file
  write_csv(microTab, "neutral_mutations_WHO_A_not_B.csv")
}

safeBinomMeld = function(x, y, z, w) {
  if (y == 0 || w == 0) {
    output = list(estimate = NA, conf.int = c(NA, NA))
  } else {
    output = binomMeld.test(x1 = x, n1 = y, x2 = z, n2 = w, parmtype = "ratio", alternative = "greater")
  }
  output
}

oldCatalogueStats = function(TabT) {
  TabT %<>%
    mutate(      LRP = pmap(list(present_S, present_S + absent_S, present_R, present_R + absent_R), safeBinomMeld)) %>%
    mutate( LRP_SOLO = pmap(list(   SOLO_S,    SOLO_S + absent_S,    SOLO_R,    SOLO_R + absent_R), safeBinomMeld)) %>%
    mutate(      LRN = pmap(list( absent_S, present_S + absent_S,  absent_R, present_R + absent_R), safeBinomMeld)) %>%
    mutate( LRN_SOLO = pmap(list( absent_S,    SOLO_S + absent_S,  absent_R,    SOLO_R + absent_R), safeBinomMeld))
}

postprocessGradedStats = function(fast = FALSE, correct_all = TRUE, skipBDQfromSA = FALSE) {
  matchTab = read_csv("new_variant_matched_to_old.csv", guess_max = Inf, show_col_types = FALSE)
  for (name in unique(PHENO_GROUPS$group)) {
    if (!fast || name == "MAIN") {
      if (name == "MAIN") {
        Tab0 = read_csv("Graded_Stats_WHO.csv" , guess_max = Inf, show_col_types = FALSE)
        Tab1 = read_csv("Graded_Stats_MAIN.csv", guess_max = Inf, show_col_types = FALSE)
        Tab01 = full_join(Tab0, Tab1, by = c("drug", "variant")) %>%
          mutate(useWHO = or(useWHO.x, useWHO.y), neutral = or(neutral.x, neutral.y)) %>%
          select(-useWHO.x, -useWHO.y, -neutral.x, -neutral.y)
        Tab01 %<>%
          mutate_at("useWHO", ~{na_if(., FALSE)}) %>%
          group_by(drug, variant) %>%
          mutate(useWHO = any(useWHO)) %>%
          ungroup() %>%
          mutate(useWHO = useWHO | neutral | (!is.na(Initial.x) & Initial.x != 3 & Initial.y == 3)) %>%
          mutate_at("useWHO",
                    ~{ifelse(!is.na(Initial.y) & (Initial.x == Initial.y | Initial.y != 3) & is.na(.), FALSE, .)}) %>%
          mutate(special_case = is.na(useWHO)) %>%
          mutate_at("useWHO", ~{replace_na(., FALSE)}) %>%
          select(drug, variant, useWHO, special_case) %>%
          group_by(drug, variant) %>%
          slice(1) %>%
          ungroup()
        useWHO = Tab01 %>% 
          dplyr::filter(useWHO)
        useALL = Tab01 %>%
          dplyr::filter(!useWHO)
        finalTab = bind_rows(Tab0 %>% select(-useWHO) %>% inner_join(useWHO, by = c("drug", "variant")), 
                             Tab1 %>% select(-useWHO) %>% inner_join(useALL, by = c("drug", "variant")))
        finalTab %<>%
          mutate(datasets = ifelse(useWHO, "WHO", "ALL")) %>%
          mutate(downgrade = (special_case & Final == 1 & datasets == "ALL")) %>%
          mutate_at("Additional grading criteria_v1", ~{ifelse(downgrade, "Downgraded to interim based on WHO dataset", .)}) %>%
          mutate_at("Final"                         , ~{ifelse(downgrade, 2,                                            .)}) %>%
          select(-special_case, -downgrade) %>%
          mutate_at("Additional grading criteria_v1", ~{ifelse(Final %in% c(1, 5) & datasets == "ALL" & is.na(.), ALL_ONLY, .)}) %>%
          mutate_at("Final",                       ~{ifelse(`Additional grading criteria_v1` == ALL_ONLY, 2, .)})
      } else {
        finalTab = read_csv(paste0("Graded_Stats_", name, ".csv"), guess_max = Inf, show_col_types = FALSE)
      }
      finalTab %<>%
        mutate_at("Additional grading criteria_v1", ~{paste(., ifelse(algorithm_pass == 2 & Final %in% c(1,5), "Algorithm pass 2", ""))}) %>%
        mutate_at("Final",                       ~{ifelse(Final == 1 & algorithm_pass == 2, 2, .)}) %>%
        mutate(Initial_Confidence_Grading = GRADES[Initial], Final_Confidence_Grading = GRADES[Final]) %>%
        select(-Initial, -Final)
      finalTab %<>%
        rename("variant_category_v2" = variant) %>%
        left_join(matchTab, by = c("drug", "variant_category_v2")) %>%
        mutate(Present_in_Catalogue_v1 = !is.na(variant_category_v1)) %>%
        rename("variant" = variant_category_v2, "variant_v1_nomenclature" = variant_category_v1)
      write_csv(finalTab, paste0(paste("Final_graded_algorithm_catalogue", name, Sys.Date(), "Leonid", 
                                       "Correct", ifelse(correct_all, "All", "SOLO"), ifelse(skipBDQfromSA, "withoutSA", ""), sep = "_"), ".csv"))
    }
  }
}

EXTRA_LOF = FALSE # Set to TRUE if you want to generate information for the other two types of pooled indels as well!
if (EXTRA_LOF) {
  POOLED_EFFECTS[["inframe"]] = paste0("inframe_", c("insertion", "deletion"))
  POOLED_EFFECTS[["LoF_ALL"]] = c(POOLED_EFFECTS[["LoF"]], POOLED_EFFECTS[["inframe"]])
}

GROUP_123_DRUGS  = DRUG_LIST[match(c("BDQ", "INH", "PZA")              , SHORT_NAMES)]

if (stage == 2 & relax == FALSE & version == 2 & geno == "ALL") {
  miniSubset = curSubset %>%
    filter(selected & phenotype == "R" & drug == "CFZ") %>%
    distinct(sample_id, .keep_all = TRUE)
  write_csv(miniSubset, "Group2_75_RIFALL_CFZ_TP_samples.csv")
}

geneMap = fullDataset %>% select(gene, drug) %>% distinct()

## Sanity check: most of the isolates containing a URM should have a resistant phenotype
## miniTab = table(masterTab %>% slice(1) %>% pull(phenotype), masterTab %>% slice(1) %>% pull(urm))
## print("Displaying the table of URM presence/absence by phenotype; one unit is a sample-drug combination")
## print(miniTab)

## From the NewCrossCheck function:
# resultsJoint %<>%
#   testAndRemove("Sens", "Sensitivity")
# resultsJoint %<>%
#   testAndRemove("sens_all_lb", "Sens_lb")
# resultsJoint %<>%
#   testAndRemove("sens_all_ub", "Sens_ub")
# resultsJoint %<>%
#   testAndRemove("Spec", "Specificity")
# resultsJoint %<>%
#   testAndRemove("spec_all_lb", "Spec_lb")
# resultsJoint %<>%
#   testAndRemove("spec_all_ub", "Spec_ub")
