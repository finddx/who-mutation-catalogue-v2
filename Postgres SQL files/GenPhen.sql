-- Country table
-- Filed using pycountry module from python
-- The primary key is the iso integer value
CREATE TABLE public."country" (
  "country_id" int,
  "two_letters_code" char(2),
  "three_letters_code" char(3),
  "country_usual_name" varchar,
  "country_official_name" varchar
);

CREATE SCHEMA "genphensql";

SET search_path TO "genphensql", biosql, public;

-- The sample table accomodates non compulsory information from NCBI BioSample
-- i.e. BioSampleId for use of the entrez API and the synonym SRAName
-- Some information from BioSample (sampling date, country of origin) is incorporated
CREATE TABLE "sample" (
  "sample_id" bigserial,
  "biosample_id" int,
  "sample_name" varchar NOT NULL,
  "sra_name" varchar,
  "ncbi_taxon_id" int NOT NULL,
  "submission_date" date,
  "sampling_date" daterange,
  "country_id" int,
  "additional_geographical_information" varchar,
  "latitude" varchar,
  "longitude" varchar,
  "isolation_source" varchar,
  "status" varchar,
  "last_status_change" timestamp,
  "patiend_id" bigint
);

-- Staged Sample table because AWS Glue PySpark does not handle daterange
-- We need to perform some UPSERT as well (when updating entries originating from NCBI BioSample)
CREATE TABLE "staged_sample" (
  "sample_name" varchar NOT NULL,
  "ncbi_taxon_id" int,
  "submission_date" varchar,
  "sampling_date" varchar,
  "country_id" int,
  "additional_geographical_information" varchar,
  "latitude" varchar,
  "longitude" varchar,
  "isolation_source" varchar,
  "patiend_id" bigint,
  "library_name" varchar,
  "prefix" varchar,
  "library_preparation_strategy" varchar,
  "sequencing_platform" varchar,
  "library_layout" varchar,
  "aws_bucket_region" varchar
);

-- Incoporating XML BioProject information as is from NCBI entrez API
CREATE TABLE "bioproject" (
  "bioproject_id" bigint,
  "ncbi_xml_value" xml
);

-- For BioProject, some informations are parsed from the XML and inserted in a secondary table
CREATE TABLE "dataset" (
  "dataset_id" bigserial,
  "bioproject_id" bigint,
  "dataset_name" varchar NOT NULL,
  "dataset_origin" varchar,
  "dataset_owner" varchar,
  "contact_email" varchar,
  "dataset_title" varchar,
  "description" varchar,
  "submission_date" date
);

-- Staging datasets for submission based on pre existing data in NCBI/ENA
CREATE TABLE "staged_dataset" (
  "dataset_name" varchar NOT NULL,
  "dataset_origin" varchar,
  "dataset_owner" varchar,
  "contact_email" varchar,
  "dataset_title" varchar,
  "description" varchar
);

-- Many Dataset <-> Many Sample
CREATE TABLE "dataset_to_sample" (
  "dataset_id" bigint,
  "sample_id" bigint
);

CREATE TABLE "contributor" (
  "contributor_id" serial,
  "contributor_name" varchar NOT NULL,
  "contributor_affiliation" varchar,
  "contributor_email" varchar
);

-- Many Dataset <-> Many Contributor
CREATE TABLE "dataset_to_contributor" (
  "dataset_id" bigint,
  "contributor_id" bigint
);

-- One Patient <-> Many Sample
CREATE TABLE "patient" (
  "patient_id" bigserial,
  "gender" char(1),
  "age_at_sampling" int,
  "disease" varchar,
  "new_tuberculosis_case" boolean,
  "previous_treatment_category" varchar,
  "treatment_regimen" varchar,
  "treatment_duration" int,
  "treatment_outcome" varchar,
  "hiv_positive" boolean
);

-- Information about the raw sequencing data. Accomodates SRA database from NCBI
-- Sequencing data can be private (data location is S3) or public (data location = NCBI)
CREATE TABLE "sequencing_data" (
  "sequencing_data_id" bigserial,
  "sample_id" int NOT NULL,
  "library_name" varchar NOT NULL,
  "data_location" varchar NOT NULL,
  "library_preparation_strategy" varchar,
  "dna_source" varchar,
  "dna_selection" varchar,
  "sequencing_platform" varchar,
  "sequencing_machine" varchar,
  "library_layout" varchar,
  "file_path" varchar,
  "assay" varchar
);

CREATE TABLE "sequencing_data_hash" (
  "sequencing_data_id" bigint NOT NULL,
  "algorithm" varchar,
  "value" varchar NOT NULL
);

-- Create another table to store sample label synonym
-- Accomodates data from BioSamples
CREATE TABLE "additional_sample_name" (
  "sample_id" int NOT NULL,
  "db" varchar NOT NULL,
  "db_label" varchar NOT NULL,
  "sample_name_synonym" varchar NOT NULL
);

-- Storing the unique coordinates necessary for description of a variant,
-- Reference name, position, reference nucleotide, alternative nucleotide
CREATE TABLE "variant" (
  "variant_id" bigserial,
  "chromosome" varchar NOT NULL,
  "position" int NOT NULL,
  "reference_nucleotide" varchar NOT NULL,
  "alternative_nucleotide" varchar NOT NULL
);

-- Storing HGVS annotations with respect to Reference Sequence Ids from BioSQL
-- References can be genes (dbxref IDs to NCBI Gene db), when annotation refers to nucleotides
-- References can be proteins (dbxref IDS to NCBI Protein db) when annotation refers to amino acid
CREATE TABLE "annotation" (
  "annotation_id" bigserial,
  "reference_db_crossref_id" int NOT NULL,
  "hgvs_value" varchar NOT NULL,
  "predicted_effect" varchar NOT NULL,
  "distance_to_reference" int
);

-- Many Variant <-> Many Annotation
-- One variant can have many annotations: annotations on different nearby genes, annotation on gene sequence or protein sequence for the same locus
-- One annotation can have many variant: redundancy of the genetic code, or more complex variant leading to same change
CREATE TABLE "variant_to_annotation" (
  "variant_id" int NOT NULL,
  "annotation_id" int NOT NULL
);

-- Staging table for both link between variant id and new annotation values
-- All new incoming data is copied here before final insertion
-- Only new annotations and new links between variants and new anotations are inserted
-- For VariantToAnnotation insertion, JOINS are performed to insert the Ids of each componenent
CREATE TABLE "staged_variant_to_annotation" (
  "variant_id" int NOT NULL,
  "locus_tag_name" varchar NOT NULL,
  "hgvs_value" varchar NOT NULL,
  "predicted_effect" varchar NOT NULL,
  "type" varchar NOT NULL,
  "distance_to_reference" int
  );

-- Storing genotype results for each samples. 
-- Only stores the Id of the Variant instead of full coordinates
-- Many Genotypes <-> One Variant
-- Many Genotypes <-> One Sample
CREATE TABLE "genotype" (
  "genotype_id" bigserial,
  "sample_id" int NOT NULL,
  "variant_id" int NOT NULL,
  "genotyper" varchar NOT NULL,
  "quality" float NOT NULL,
  "reference_ad" int NOT NULL,
  "alternative_ad" int NOT NULL,
  "total_dp" int NOT NULL,
  "genotype_value" varchar NOT NULL
);

-- Storing the results of the PySpark logic for the variant categories
CREATE TABLE "tiered_variant_categories" (
  "variant_id" int NOT NULL,
  "gene_db_crossref_id" int NOT NULL,
  "hgvs_value" varchar NOT NULL,
  "predicted_effect" varchar NOT NULL,
  "distance_to_reference" int
);

-- Storing PCA results
CREATE TABLE "pca_dimension_results" (
  "sample_id" int NOT NULL,
  "dimension" int NOT NULL,
  "value" float NOT NULL
);

-- Staging table for Genotypes
-- Incomming data is copied here then inserted after JOINING with "variant"
CREATE TABLE "staged_genotype" (
  "sample_name" varchar NOT NULL,
  "chromosome" varchar NOT NULL,
  "position" int NOT NULL,
  "variant_id" int, 
  "reference_nucleotide" varchar NOT NULL,
  "alternative_nucleotide" varchar NOT NULL,
  "genotyper" varchar NOT NULL,
  "quality" float NOT NULL,
  "reference_ad" int NOT NULL,
  "alternative_ad" int NOT NULL,
  "total_dp" int NOT NULL,
  "genotype_value" varchar NOT NULL
);

-- Stores the sequencing depth of each locus for each sample
CREATE TABLE "locus_sequencing_stats" (
  "sample_id" int NOT NULL,
  "gene_db_crossref_id" int NOT NULL,
  "mean_depth" float NOT NULL,
  "coverage_10x" float NOT NULL,
  "coverage_15x" float NOT NULL,
  "coverage_20x" float NOT NULL,
  "coverage_30x" float NOT NULL
);

-- High level statistics of read mapping onto the reference genome per sample
CREATE TABLE "summary_sequencing_stats" (
  "sample_id" int NOT NULL,
  "median_depth" float NOT NULL,
  "coverage_10x" float NOT NULL,
  "coverage_15x" float NOT NULL,
  "coverage_20x" float NOT NULL,
  "coverage_30x" float NOT NULL,
  "raw_total_sequences" bigint NOT NULL,
  "filtered_sequences" bigint NOT NULL,
  "sequences" bigint NOT NULL,
  "is_sorted" bigint NOT NULL,
  "first_fragments" bigint NOT NULL,
  "last_fragments" bigint NOT NULL,
  "reads_mapped" bigint NOT NULL,
  "reads_mapped_and_paired" bigint NOT NULL,
  "reads_unmapped" bigint NOT NULL,
  "reads_properly_paired" bigint NOT NULL,
  "reads_paired" bigint NOT NULL,
  "reads_duplicated" bigint NOT NULL,
  "reads_mq_0" bigint NOT NULL,
  "reads_qc_failed" bigint NOT NULL,
  "non_primary_alignments" bigint NOT NULL,
  "total_length" bigint NOT NULL,
  "total_first_fragment_length" bigint NOT NULL,
  "total_last_fragment_length" bigint NOT NULL,
  "bases_mapped" bigint NOT NULL,
  "bases_mapped_cigar" bigint NOT NULL,
  "bases_trimmed" bigint NOT NULL,
  "bases_duplicated" bigint NOT NULL,
  "mismatches" bigint NOT NULL,
  "error_rate" float NOT NULL,
  "average_length" int NOT NULL,
  "average_first_fragment_length" int NOT NULL,
  "average_last_fragment_length" int NOT NULL,
  "maximum_length" int NOT NULL,
  "maximum_first_fragment_length" int NOT NULL,
  "maximum_last_fragment_length" int NOT NULL,
  "average_quality" float NOT NULL,
  "insert_size_average" float NOT NULL,
  "insert_size_standard_deviation" float NOT NULL,
  "inward_oriented_pairs" int NOT NULL,
  "outward_oriented_pairs" int NOT NULL,
  "pairs_with_other_orientation" int NOT NULL,
  "pairs_on_different_chromosomes" int NOT NULL,
  "percentage_of_properly_paired_reads" float NOT NULL
);

-- Stores taxonomy information extract from kraken analysis of raw reads
CREATE TABLE "taxonomy_stats" (
  "sample_id" int NOT NULL,
  "ncbi_taxon_id" int NOT NULL,
  "value" float NOT NULL
);

-- Target of the targeted Next Generation Sequencing kits
CREATE TABLE "amplicon_target" (
  "amplicon_target_id" serial,
  "amplicon_assay_name" varchar,
  "chromosome" varchar NOT NULL,
  "start" int NOT NULL,
  "end" int  NOT NULL,
  "gene_db_crossref_id" int NOT NULL
);

CREATE TABLE "growth_medium" (
  "medium_id" serial,
  "medium_name" varchar NOT NULL UNIQUE
);

CREATE TABLE "phenotypic_drug_susceptibility_assessment_method" (
  "method_id" serial,
  "method_name" varchar NOT NULL UNIQUE
);

CREATE TABLE "drug" (
  "drug_id" serial,
  "drug_name" varchar NOT NULL UNIQUE
);

CREATE TABLE "drug_synonym" (
  "drug_id" bigint,
  "code" varchar,
  "drug_name_synonym" varchar NOT NULL UNIQUE
);

CREATE TABLE "microdilution_plate_concentration" (
  "plate" varchar NOT NULL,
  "drug_id" int NOT NULL,
  "concentration" float NOT NULL 
);

-- Associations between gene and drug for resistance markers
CREATE TABLE "gene_drug_resistance_association" (
  "gene_db_crossref_id" int NOT NULL,
  "drug_id" int NOT NULL,
  "tier" int NOT NULL
);

CREATE TABLE "phenotypic_drug_susceptibility_test"(
   "test_id" bigserial,
   "sample_id" int NOT NULL,
   "drug_id"  int NOT NULL,
   "medium_id" int,
   "method_id" int,
   "concentration" float,
   "test_result" char(1) NOT NULL,
   "submission_date" date
);

CREATE TABLE "minimum_inhibitory_concentration_test"(
   "test_id" bigserial,
   "sample_id" int NOT NULL,
   "drug_id"  int NOT NULL,
   "plate" varchar,
   "mic_value" numrange NOT NULL,
   "submission_date" date

);

-- Staged MIC table because AWS Glue PySpark does not handle numrange
CREATE TABLE "staged_minimum_inhibitory_concentration_test"(
   "sample_id" int NOT NULL,
   "drug_id"  int NOT NULL,
   "plate" varchar,
   "mic_value" varchar NOT NULL
);

CREATE TABLE "smear_microscopy_results"(
  "test_id" bigserial,
  "sample_id" bigint NOT NULL,
  "smear_result" varchar
);

CREATE TABLE "molecular_drug_resistance_test"(
  "test_id" bigserial,
  "sample_id" bigint NOT NULL,
  "test_name" varchar NOT NULL,
  "drug_id" int NOT NULL,
  "test_result" char(1) NOT NULL
);

CREATE TABLE "epidemiological_cut_off_value" (
    "drug_id" integer NOT NULL,
    "medium_name" varchar NOT NULL,
    "value" float NOT NULL
);

INSERT INTO "drug"("drug_name")
SELECT x."Drug"
FROM (
  VALUES 
  ('H', 'Isoniazid'),
  ('R', 'Rifampicin'),
  ('S', 'Streptomycin'),
  ('E', 'Ethambutol'),
  ('Ofx', 'Ofloxacin'),
  ('Cm', 'Capreomycin'),
  ('Am', 'Amikacin'),
  ('Km', 'Kanamycin'),
  ('Z', 'Pyrazinamide'),
  ('Lfx', 'Levofloxacin'),
  ('Mfx', 'Moxifloxacin'),
  ('Pas', 'Para-Aminosalicylic Acid'),
  ('Pto', 'Prothionamide'),
  ('Cs', 'Cycloserine'),
  ('Amx/Clv','Amoxicillin-Clavulanate'),
  ('Mb', 'Rifabutin'),
  ('Eto', 'Ethionamide'),
  ('Dld', 'Delamanid'),
  ('Bdq', 'Bedaquiline'),
  ('Ipm/Cln', 'Imipenem-Cilastatin'),
  ('Lzd', 'Linezolid'),
  ('Cfz', 'Clofazimine'),
  ('Clr', 'Clarithromycin'),
  ('Ft', 'Fluoroquinolones'),
  ('AG/CP', 'Aminoglycosides'),
  ('Gfx', 'Gatifloxacin'),
  ('Cip', 'Ciprofloxacin'),
  ('Sit', 'Sitafloxacin'),
  ('Azt', 'Azithromycin'),
  ('PA-284', 'Pretomanid')
  ) x("Synonym", "Drug")
;

INSERT INTO "drug_synonym"("drug_id", "drug_name_synonym")
SELECT "drug"."drug_id", x."Synonym"
FROM (
  VALUES 
  ('H', 'Isoniazid'),
  ('R', 'Rifampicin'),
  ('S', 'Streptomycin'),
  ('E', 'Ethambutol'),
  ('Ofx', 'Ofloxacin'),
  ('Cm', 'Capreomycin'),
  ('Am', 'Amikacin'),
  ('Km', 'Kanamycin'),
  ('Z', 'Pyrazinamide'),
  ('Lfx', 'Levofloxacin'),
  ('Mfx', 'Moxifloxacin'),
  ('Pas', 'Para-Aminosalicylic Acid'),
  ('Pto', 'Prothionamide'),
  ('Cs', 'Cycloserine'),
  ('Dcs', 'Cycloserine'),
  ('Amx/Clv','Amoxicillin-Clavulanate'),
  ('Mb', 'Rifabutin'),
  ('Eto', 'Ethionamide'),
  ('Dld', 'Delamanid'),
  ('Bdq', 'Bedaquiline'),
  ('Ipm/Cln', 'Imipenem-Cilastatin'),
  ('Lzd', 'Linezolid'),
  ('Cfz', 'Clofazimine'),
  ('Clr', 'Clarithromycin'),
  ('Ft', 'Fluoroquinolones'),
  ('AG/CP', 'Aminoglycosides'),
  ('Gfx', 'Gatifloxacin'),
  ('Cip', 'Ciprofloxacin'),
  ('Sit', 'Sitafloxacin'),
  ('Azt', 'Azithromycin'),
  ('Pa', 'Pretomanid'),
  ('PA-824', 'Pretomanid'),
  ('Mycobutin', 'Rifabutin')
  ) x("Synonym", "Drug")
INNER JOIN "drug" on "drug"."drug_name"=x."Drug";

INSERT INTO "drug_synonym"("drug_id", "code", "drug_name_synonym")
SELECT "drug"."drug_id",
      'three_letter_code',
      x."Synonym"
FROM (
  VALUES
    ('INH', 'Isoniazid'),
    ('RIF', 'Rifampicin'),
    ('STR', 'Streptomycin'),
    ('STM', 'Streptomycin'),
    ('EMB', 'Ethambutol'),
    ('PZA', 'Pyrazinamide'),
    ('LFX', 'Levofloxacin'),
    ('LVX', 'Levofloxacin'),
    ('LEV', 'Levofloxacin'),
    ('LEVO', 'Levofloxacin'),
    ('MFX', 'Moxifloxacin'),
    ('MXF', 'Moxifloxacin'),
    ('MOXI', 'Moxifloxacin'),
    ('MOX', 'Moxifloxacin'),    
    ('OFX', 'Ofloxacin'),
    ('OFL', 'Ofloxacin'),
    ('CAP', 'Capreomycin'),
    ('AMK', 'Amikacin'),
    ('AMI', 'Amikacin'),
    ('KAN', 'Kanamycin'),
    ('ETH', 'Ethionamide'),
    ('PAS', 'Para-Aminosalicylic Acid'),
    ('PTO', 'Prothionamide'),
    ('DLM', 'Delamanid'),
    ('BDQ', 'Bedaquiline'),
    ('LZD', 'Linezolid'),
    ('CFZ', 'Clofazimine'),
    ('CLR', 'Clarithromycin'),
    ('GFX', 'Gatifloxacin'),
    ('CIP', 'Ciprofloxacin'),
    ('SIT', 'Sitafloxacin'),
    ('STX', 'Sitafloxacin'),
    ('AZM', 'Azithromycine'),
    ('GEN', 'Gentamicin'),
    ('MPM', 'Meropenem'),
    ('IMI', 'Imipenem'),
    ('CYC', 'Cycloserine'),
    ('DCS', 'Cycloserine'),
    ('RFB', 'Rifabutin'),
    ('PMD', 'Pretomanid')
) x("Synonym", "Drug")
INNER JOIN "drug" on "drug"."drug_name"=x."Drug";

INSERT INTO "microdilution_plate_concentration"("plate", "drug_id", "concentration")
SELECT "Plate",
      "drug_synonym"."drug_id",
      "Concentration"
FROM (
  VALUES
    ('UKMYC5', 'AMI', 8),
    ('UKMYC5', 'AMI', 4),
    ('UKMYC5', 'AMI', 2),
    ('UKMYC5', 'AMI', 1),
    ('UKMYC5', 'AMI', 0.5),
    ('UKMYC5', 'AMI', 0.25),
    ('UKMYC5', 'BDQ', 2),
    ('UKMYC5', 'BDQ', 1),
    ('UKMYC5', 'BDQ', 0.5),
    ('UKMYC5', 'BDQ', 0.25),
    ('UKMYC5', 'BDQ', 0.12),
    ('UKMYC5', 'BDQ', 0.06),
    ('UKMYC5', 'BDQ', 0.03),
    ('UKMYC5', 'BDQ', 0.015),
    ('UKMYC5', 'CFZ', 4),
    ('UKMYC5', 'CFZ', 2),
    ('UKMYC5', 'CFZ', 1),
    ('UKMYC5', 'CFZ', 0.5),
    ('UKMYC5', 'CFZ', 0.25),
    ('UKMYC5', 'CFZ', 0.12),
    ('UKMYC5', 'CFZ', 0.06),
    ('UKMYC5', 'DLM', 1),
    ('UKMYC5', 'DLM', 0.5),
    ('UKMYC5', 'DLM', 0.25),
    ('UKMYC5', 'DLM', 0.12),
    ('UKMYC5', 'DLM', 0.06),
    ('UKMYC5', 'DLM', 0.03),
    ('UKMYC5', 'DLM', 0.015),
    ('UKMYC5', 'EMB', 8),
    ('UKMYC5', 'EMB', 4),
    ('UKMYC5', 'EMB', 2),
    ('UKMYC5', 'EMB', 1),
    ('UKMYC5', 'EMB', 0.5),
    ('UKMYC5', 'EMB', 0.25),
    ('UKMYC5', 'EMB', 0.06),
    ('UKMYC5', 'EMB', 0.12),
    ('UKMYC5', 'ETH', 8),
    ('UKMYC5', 'ETH', 4),
    ('UKMYC5', 'ETH', 2),
    ('UKMYC5', 'ETH', 1),
    ('UKMYC5', 'ETH', 0.5),
    ('UKMYC5', 'ETH', 0.25),
    ('UKMYC5', 'INH', 1.6),
    ('UKMYC5', 'INH', 0.8),
    ('UKMYC5', 'INH', 0.4),
    ('UKMYC5', 'INH', 0.2),
    ('UKMYC5', 'INH', 0.1),
    ('UKMYC5', 'INH', 0.05),
    ('UKMYC5', 'INH', 0.025),
    ('UKMYC5', 'KAN', 16),
    ('UKMYC5', 'KAN', 8),
    ('UKMYC5', 'KAN', 4),
    ('UKMYC5', 'KAN', 2),
    ('UKMYC5', 'KAN', 1),
    ('UKMYC5', 'LEV', 8),
    ('UKMYC5', 'LEV', 4),
    ('UKMYC5', 'LEV', 2),
    ('UKMYC5', 'LEV', 1),
    ('UKMYC5', 'LEV', 0.5),
    ('UKMYC5', 'LEV', 0.25),
    ('UKMYC5', 'LEV', 0.12),
    ('UKMYC5', 'LZD', 2),
    ('UKMYC5', 'LZD', 1),
    ('UKMYC5', 'LZD', 0.5),
    ('UKMYC5', 'LZD', 0.25),
    ('UKMYC5', 'LZD', 0.12),
    ('UKMYC5', 'LZD', 0.06),
    ('UKMYC5', 'LZD', 0.03),
    ('UKMYC5', 'MXF', 4),
    ('UKMYC5', 'MXF', 2),
    ('UKMYC5', 'MXF', 1),
    ('UKMYC5', 'MXF', 0.5),
    ('UKMYC5', 'MXF', 0.25),
    ('UKMYC5', 'MXF', 0.12),
    ('UKMYC5', 'MXF', 0.06),
    ('UKMYC5', 'PAS', 4),
    ('UKMYC5', 'PAS', 2),
    ('UKMYC5', 'PAS', 1),
    ('UKMYC5', 'PAS', 0.5),
    ('UKMYC5', 'PAS', 0.25),
    ('UKMYC5', 'PAS', 0.12),
    ('UKMYC5', 'RFB', 2),
    ('UKMYC5', 'RFB', 1),
    ('UKMYC5', 'RFB', 0.5),
    ('UKMYC5', 'RFB', 0.25),
    ('UKMYC5', 'RFB', 0.12),
    ('UKMYC5', 'RFB', 0.06),
    ('UKMYC5', 'RIF', 4),
    ('UKMYC5', 'RIF', 2),
    ('UKMYC5', 'RIF', 1),
    ('UKMYC5', 'RIF', 0.5),
    ('UKMYC5', 'RIF', 0.25),
    ('UKMYC5', 'RIF', 0.12),
    ('UKMYC5', 'RIF', 0.06),
    ('UKMYC6', 'AMI', 16),
    ('UKMYC6', 'AMI', 8),
    ('UKMYC6', 'AMI', 4),
    ('UKMYC6', 'AMI', 2),
    ('UKMYC6', 'AMI', 1),
    ('UKMYC6', 'AMI', 0.5),
    ('UKMYC6', 'AMI', 0.25),
    ('UKMYC6', 'BDQ', 1),
    ('UKMYC6', 'BDQ', 0.5),
    ('UKMYC6', 'BDQ', 0.25),
    ('UKMYC6', 'BDQ', 0.12),
    ('UKMYC6', 'BDQ', 0.06),
    ('UKMYC6', 'BDQ', 0.03),
    ('UKMYC6', 'BDQ', 0.015),
    ('UKMYC6', 'BDQ', 0.008),
    ('UKMYC6', 'CFZ', 2),
    ('UKMYC6', 'CFZ', 1),
    ('UKMYC6', 'CFZ', 0.5),
    ('UKMYC6', 'CFZ', 0.25),
    ('UKMYC6', 'CFZ', 0.12),
    ('UKMYC6', 'CFZ', 0.06),
    ('UKMYC6', 'CFZ', 0.03),
    ('UKMYC6', 'DLM', 0.5),
    ('UKMYC6', 'DLM', 0.25),
    ('UKMYC6', 'DLM', 0.12),
    ('UKMYC6', 'DLM', 0.06),
    ('UKMYC6', 'DLM', 0.03),
    ('UKMYC6', 'DLM', 0.015),
    ('UKMYC6', 'DLM', 0.008),
    ('UKMYC6', 'EMB', 32),
    ('UKMYC6', 'EMB', 16),
    ('UKMYC6', 'EMB', 8),
    ('UKMYC6', 'EMB', 4),
    ('UKMYC6', 'EMB', 2),
    ('UKMYC6', 'EMB', 1),
    ('UKMYC6', 'EMB', 0.5),
    ('UKMYC6', 'EMB', 0.25),
    ('UKMYC6', 'ETH', 8),
    ('UKMYC6', 'ETH', 4),
    ('UKMYC6', 'ETH', 2),
    ('UKMYC6', 'ETH', 1),
    ('UKMYC6', 'ETH', 0.5),
    ('UKMYC6', 'ETH', 0.25),
    ('UKMYC6', 'INH', 12.8),
    ('UKMYC6', 'INH', 6.4),
    ('UKMYC6', 'INH', 3.2),
    ('UKMYC6', 'INH', 1.6),
    ('UKMYC6', 'INH', 0.8),
    ('UKMYC6', 'INH', 0.4),
    ('UKMYC6', 'INH', 0.2),
    ('UKMYC6', 'INH', 0.1),
    ('UKMYC6', 'INH', 0.05),
    ('UKMYC6', 'INH', 0.025),
    ('UKMYC6', 'KAN', 16),
    ('UKMYC6', 'KAN', 8),
    ('UKMYC6', 'KAN', 4),
    ('UKMYC6', 'KAN', 2),
    ('UKMYC6', 'KAN', 1),
    ('UKMYC6', 'LEV', 8),
    ('UKMYC6', 'LEV', 4),
    ('UKMYC6', 'LEV', 2),
    ('UKMYC6', 'LEV', 1),
    ('UKMYC6', 'LEV', 0.5),
    ('UKMYC6', 'LEV', 0.25),
    ('UKMYC6', 'LEV', 0.12),
    ('UKMYC6', 'LZD', 4),
    ('UKMYC6', 'LZD', 2),
    ('UKMYC6', 'LZD', 1),
    ('UKMYC6', 'LZD', 0.5),
    ('UKMYC6', 'LZD', 0.25),
    ('UKMYC6', 'LZD', 0.12),
    ('UKMYC6', 'LZD', 0.06),
    ('UKMYC6', 'MXF', 4),
    ('UKMYC6', 'MXF', 2),
    ('UKMYC6', 'MXF', 1),
    ('UKMYC6', 'MXF', 0.5),
    ('UKMYC6', 'MXF', 0.25),
    ('UKMYC6', 'MXF', 0.12),
    ('UKMYC6', 'MXF', 0.06),
    ('UKMYC6', 'RIF', 8),
    ('UKMYC6', 'RIF', 4),
    ('UKMYC6', 'RIF', 2),
    ('UKMYC6', 'RIF', 1),
    ('UKMYC6', 'RIF', 0.5),
    ('UKMYC6', 'RIF', 0.25),
    ('UKMYC6', 'RIF', 0.12),
    ('UKMYC6', 'RIF', 0.06),
    ('UKMYC6', 'RIF', 0.03),
    ('UKMYC6', 'RFB', 2),
    ('UKMYC6', 'RFB', 1),
    ('UKMYC6', 'RFB', 0.5),
    ('UKMYC6', 'RFB', 0.25),
    ('UKMYC6', 'RFB', 0.12),
    ('UKMYC6', 'RFB', 0.06),
    ('MYCOTB', 'AMI', 16),
    ('MYCOTB', 'AMI', 8),
    ('MYCOTB', 'AMI', 4),
    ('MYCOTB', 'AMI', 2),
    ('MYCOTB', 'AMI', 1),
    ('MYCOTB', 'AMI', 0.5),
    ('MYCOTB', 'AMI', 0.25),
    ('MYCOTB', 'AMI', 0.12),
    ('MYCOTB', 'CYC', 256),
    ('MYCOTB', 'CYC', 128),
    ('MYCOTB', 'CYC', 64),
    ('MYCOTB', 'CYC', 32),
    ('MYCOTB', 'CYC', 16),
    ('MYCOTB', 'CYC', 8),
    ('MYCOTB', 'CYC', 4),
    ('MYCOTB', 'CYC', 2),
    ('MYCOTB', 'EMB', 32),
    ('MYCOTB', 'EMB', 16),
    ('MYCOTB', 'EMB', 8),
    ('MYCOTB', 'EMB', 4),
    ('MYCOTB', 'EMB', 2),
    ('MYCOTB', 'EMB', 1),
    ('MYCOTB', 'EMB', 0.5),
    ('MYCOTB', 'ETH', 40),
    ('MYCOTB', 'ETH', 20),
    ('MYCOTB', 'ETH', 10),
    ('MYCOTB', 'ETH', 5),
    ('MYCOTB', 'ETH', 2.5),
    ('MYCOTB', 'ETH', 1.2),
    ('MYCOTB', 'ETH', 0.6),
    ('MYCOTB', 'ETH', 0.3),
    ('MYCOTB', 'INH', 4),
    ('MYCOTB', 'INH', 2),
    ('MYCOTB', 'INH', 1),
    ('MYCOTB', 'INH', 0.5),
    ('MYCOTB', 'INH', 0.25),
    ('MYCOTB', 'INH', 0.12),
    ('MYCOTB', 'INH', 0.06),
    ('MYCOTB', 'INH', 0.03),
    ('MYCOTB', 'KAN', 40),
    ('MYCOTB', 'KAN', 20),
    ('MYCOTB', 'KAN', 10),
    ('MYCOTB', 'KAN', 5),
    ('MYCOTB', 'KAN', 2.5),
    ('MYCOTB', 'KAN', 1.2),
    ('MYCOTB', 'KAN', 0.6),
    ('MYCOTB', 'MXF', 8),
    ('MYCOTB', 'MXF', 4),
    ('MYCOTB', 'MXF', 2),
    ('MYCOTB', 'MXF', 1),
    ('MYCOTB', 'MXF', 0.5),
    ('MYCOTB', 'MXF', 0.25),
    ('MYCOTB', 'MXF', 0.12),
    ('MYCOTB', 'MXF', 0.06),
    ('MYCOTB', 'OFL', 32),
    ('MYCOTB', 'OFL', 16),
    ('MYCOTB', 'OFL', 8),
    ('MYCOTB', 'OFL', 4),
    ('MYCOTB', 'OFL', 2),
    ('MYCOTB', 'OFL', 1),
    ('MYCOTB', 'OFL', 0.5),
    ('MYCOTB', 'OFL', 0.25),
    ('MYCOTB', 'PAS', 64),
    ('MYCOTB', 'PAS', 32),
    ('MYCOTB', 'PAS', 16),
    ('MYCOTB', 'PAS', 8),
    ('MYCOTB', 'PAS', 4),
    ('MYCOTB', 'PAS', 2),
    ('MYCOTB', 'PAS', 1),
    ('MYCOTB', 'PAS', 0.5),
    ('MYCOTB', 'RFB', 16),
    ('MYCOTB', 'RFB', 8),
    ('MYCOTB', 'RFB', 4),
    ('MYCOTB', 'RFB', 2),
    ('MYCOTB', 'RFB', 1),
    ('MYCOTB', 'RFB', 0.5),
    ('MYCOTB', 'RFB', 0.25),
    ('MYCOTB', 'RFB', 0.12),
    ('MYCOTB', 'RIF', 16),
    ('MYCOTB', 'RIF', 8),
    ('MYCOTB', 'RIF', 4),
    ('MYCOTB', 'RIF', 2),
    ('MYCOTB', 'RIF', 1),
    ('MYCOTB', 'RIF', 0.5),
    ('MYCOTB', 'RIF', 0.25),
    ('MYCOTB', 'RIF', 0.12),
    ('MYCOTB', 'STR', 32),
    ('MYCOTB', 'STR', 16),
    ('MYCOTB', 'STR', 8),
    ('MYCOTB', 'STR', 4),
    ('MYCOTB', 'STR', 2),
    ('MYCOTB', 'STR', 1),
    ('MYCOTB', 'STR', 0.5),
    ('MYCOTB', 'STR', 0.25)
) x("Plate", "Code", "Concentration")
INNER JOIN "drug_synonym" on "drug_synonym"."drug_name_synonym"=x."Code";

INSERT INTO "growth_medium"("medium_name")
VALUES 
  ('MGIT'),
  ('BACTEC460'),
  ('LJ'),
  ('Agar'),
  ('Middlebrook7H9'),
  ('Middlebrook7H10'),
  ('Middlebrook7H11'),
  ('Waynes'),
  ('Marks Biphasic'),
  ('MODS')
;

INSERT INTO "phenotypic_drug_susceptibility_assessment_method"("method_name")
VALUES
  ('Resistance Ratio'),
  ('Proportions'),
  ('Direct'),
  ('Nitrate reductase assay'),
  ('WHO')
;

INSERT INTO "epidemiological_cut_off_value"("drug_id", "medium_name", "value")
SELECT "drug"."drug_id", y."Plate", x."Value"
FROM (
    VALUES 
    ('Isoniazid', 0.1),
    ('Rifampicin', 0.5),
    ('Ethambutol', 4),
    ('Moxifloxacin', 1),
    ('Levofloxacin', 1),
    ('Kanamycin', 4),
    ('Amikacin', 1),
    ('Rifabutin', 0.12),
    ('Clofazimine', 0.25),
    ('Linezolid', 1),
    ('Delamanid', 0.12),
    ('Bedaquiline', 0.25),
    ('Ethionamide', 4)
 ) x ("Name", "Value") 
  INNER JOIN "drug" ON "drug"."drug_name"=x."Name", 
  (
    VALUES 
    ('UKMYC5'),
    ('UKMYC6')
  ) y ("Plate"); 

CREATE TABLE "phenotypic_drug_susceptiblity_test_who_category" (
    "drug_id" integer NOT NULL,
    "medium_id" integer NOT NULL,
    "concentration" float NOT NULL,
    "category" varchar NOT NULL
);

CREATE TABLE "promoter_distance" (
  "gene_db_crossref_id" integer NOT NULL,
  "region_start" integer NOT NULL,
  "region_end" integer NOT NULL
);