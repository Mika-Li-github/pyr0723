library(biomaRt)
library(BiocManager)
# listMarts()
#browseVignettes("biomaRt")

# ENSEMBL_MART_ENSEMBL <- useMart("ENSEMBL_MART_ENSEMBL")
# listDatasets(mart = ENSEMBL_MART_ENSEMBL)
# searchDatasets(mart = ENSEMBL_MART_ENSEMBL, pattern = "GRCh38.p14")
#80 hsapiens_gene_ensembl Human genes (GRCh38.p14) GRCh38.p14

ENSEMBL_MART_SNP <- useMart("ENSEMBL_MART_SNP")
# listDatasets(mart = ENSEMBL_MART_SNP)
# searchDatasets(mart = ENSEMBL_MART_SNP, pattern = "GRCh38.p14")
#12           hsapiens_snp         Human Short Variants (SNPs and indels excluding flagged variants) (GRCh38.p14) GRCh38.p14

# ensemblGenes <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# listAttributes(mart = ensemblGenes)
# searchAttributes(mart = ensemblGenes, pattern = "entrez")
# searchAttributes(mart = ensemblGenes, pattern = "ensembl_gene_")
#3021         ensembl_gene_id         Gene stable ID          snp
# searchAttributes(mart = ensemblGenes, pattern = "snp")
#較相關的如下
#                                                               name                             description        page
# 3021                                               ensembl_gene_id                          Gene stable ID         snp
# 3024                                         ensembl_transcript_id                    Transcript stable ID         snp
# 3030                                               chromosome_name                Chromosome/scaffold name         snp
# 3031                                                start_position                         Gene start (bp)         snp
# 3032                                                  end_position                           Gene end (bp)         snp
# 3033                                                        strand                                  Strand         snp
# 3034                                                          band                          Karyotype band         snp
# 3035                                            external_gene_name                               Gene name         snp
# 3036                                          external_gene_source                     Source of gene name         snp
# 3039                                                variation_name                            Variant name         snp
# 3042                                                        allele                         Variant alleles         snp
# 3045                                                  minor_allele                            Minor allele         snp
# 3046                                             minor_allele_freq                  Minor allele frequency         snp
# 3047                                            minor_allele_count                      Minor allele count         snp
# 3050                                         snp_chromosome_strand               Variant chromosome Strand         snp
# 3052                                              chromosome_start chromosome/scaffold position start (bp)         snp
# 3053                                                chromosome_end   Chromosome/scaffold position end (bp)         snp
# 3061                                                 peptide_shift                          Protein allele         snp
# 3062                                             synonymous_status                     Variant consequence         snp
# 3063                                            allele_string_2076             Consequence specific allele         snp


ensemblSNP <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# searchAttributes(mart = ensemblSNP, pattern = "snp")
#較相關的如下
#                               name                                     description      page
# 1                        refsnp_id                                    Variant name       snp
# 2                    refsnp_source                                  Variant source       snp
# 3        refsnp_source_description                      Variant source description       snp
# 4                         chr_name                        Chromosome/scaffold name       snp
# 5                      chrom_start         Chromosome/scaffold position start (bp)       snp
# 6                        chrom_end           Chromosome/scaffold position end (bp)       snp
# 7                     chrom_strand                                          Strand       snp
# 8                           allele                                 Variant alleles       snp
# 10                       validated                     Variant supporting evidence       snp
# 11                        allele_1                                Ancestral allele       snp
# 12                    minor_allele                              Minor allele (ALL)       snp
# 13               minor_allele_freq Global minor allele frequency (all individuals)       snp
# 14              minor_allele_count     Global minor allele count (all individuals)       snp
# 15           clinical_significance                           Clinical significance       snp
# 16                    synonym_name                                    Synonym name       snp
# 17                  synonym_source                                  Synonym source       snp
# 18      synonym_source_description                      Synonym source description       snp
# 19                 variation_names                        Associated variant names       snp
# 20                      study_name                                      Study name       snp
# 21                      study_type                                      Study type       snp
# 22              study_external_ref                        Study External Reference       snp
# 23               study_description                               Study Description       snp
# 24                     source_name                                     Source name       snp
# 25                 associated_gene                  Associated gene with phenotype       snp
# 26                  phenotype_name                                  Phenotype name       snp
# 27           phenotype_description                           Phenotype description       snp
# 28  associated_variant_risk_allele                  Associated variant risk allele       snp
# 29                         p_value                                         P value       snp
# 30                        set_name                                Variant Set Name       snp
# 31                 set_description                         Variant Set Description       snp
# 32                           title                                           Title       snp
# 33                         authors                                         Authors       snp
# 34                            year                                            Year       snp
# 35                            pmid                                       PubMed ID       snp
# 36                           pmcid                    PMC reference number (PMCID)       snp
# 37                         ucsc_id                                         UCSC ID       snp
# 38                             doi                       Digital Object Identifier       snp
# 39          ensembl_gene_stable_id                                  Gene stable ID       snp
# 40               ensembl_gene_name                                       Gene Name       snp
# 42 ensembl_transcript_chrom_strand                               Transcript strand       snp
# 43                    ensembl_type                                         Biotype       snp

# listFilters(mart = ensemblSNP)
#較為實用:
#11                snp_filter                              Filter by Variant mame (e.g. rs123, CM000001)

SNP_table <- getBM(
  mart = ensemblSNP,
  filters = "snp_filter", 
  values = c("rs3124998", "rs10995271"), 
  attributes = c(
    "refsnp_id",
    "refsnp_source",
    "refsnp_source_description",
    "chr_name",
    "chrom_start",
    "chrom_end",
    "chrom_strand",
    "allele",
    "mapweight",
    "validated",
    "allele_1",
    "minor_allele",
    "minor_allele_freq",
    "minor_allele_count",
    "clinical_significance",
    "synonym_name",
    "synonym_source",
    "synonym_source_description",
    "variation_names",
    "study_name",
    "study_type",
    "study_external_ref",
    "study_description",
    "source_name",
    "associated_gene",
    "phenotype_name",
    "phenotype_description",
    "associated_variant_risk_allele",
    "p_value",
    "set_name",
    "set_description",
    "title",
    "authors",
    "year",
    "pmid",
    "pmcid",
    "ucsc_id",
    "doi",
    "ensembl_gene_stable_id",
    "ensembl_gene_name",
    "ensembl_transcript_stable_id",
    "ensembl_transcript_chrom_strand",
    "ensembl_type",
    "consequence_type_tv",
    "consequence_allele_string",
    "ensembl_peptide_allele",
    "cdna_start",
    "cdna_end",
    "distance_to_transcript",
    "polyphen_prediction",
    "polyphen_score",
    "sift_prediction",
    "sift_score",
    "reg_feature_stable_id",
    "reg_allele_string",
    "reg_consequence_types",
    "motif_feature_stable_id",
    "motif_allele_string",
    "motif_consequence_types",
    "motif_in_informative_position",
    "motif_score_delta",
    "motif_name",
    "motif_start",
    "snp"
  )
)
SNP_table
#把要的挑出來
# "refsnp_id",
# "refsnp_source",
# "chr_name",
# "chrom_start",
# "chrom_end",
# "chrom_strand",
# "allele",
# "allele_1",
# "minor_allele",
# "minor_allele_freq",
# "minor_allele_count",
# "variation_names",
# "source_name",
# "associated_gene",
# "associated_variant_risk_allele",
# "p_value",
# "ensembl_gene_stable_id",
# "ensembl_gene_name",
# "ensembl_transcript_stable_id",
# "ensembl_transcript_chrom_strand",
# "ensembl_type",
# "consequence_type_tv",
# "consequence_allele_string",
# "ensembl_peptide_allele",
# "distance_to_transcript",
# "snp"



# 
# results <- getBM(attributes = c(
#   "ensembl_gene_name",  #entrenz id，沒有的可能要去tw biobank找，population找EAS就好
#   "allele", 
#   "minor_allele", 
#   "minor_allele_freq"
# ), 
# filters = "snp_filter", 
# values = "rs4552569", 
# mart = snp)
# results
# 

#MAF, ref allele, minor allele搞定(rs71559680)


