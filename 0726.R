#####20240726#######
library(biomaRt)
library(BiocManager)
#listEnsembl()
#genes和snps可能是我們要的
#先看genes
ensembl_g <- useEnsembl(biomart = "genes")
#searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl_g <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl_g)

k <- searchAttributes(mart = ensembl_g, pattern = "snp")
#再來看snps
ensembl_s <- useEnsembl(biomart = "snps")
#searchDatasets(mart = ensembl_s, pattern = "hsapiens")
#發現有四個，選擇hsapiens_snp
ensembl_s <- useDataset(dataset = "hsapiens_snp", mart = ensembl_s)
listAttributes(ensembl_s)
j <- searchAttributes(mart = ensembl_s, pattern = "snp")

# 找出在 df1 中但不在 df2 中的 ID
unique_to_k <- setdiff(k$name, j$name)

# 找出在 df2 中但不在 df1 中的 ID
unique_to_j <- setdiff(j$name, k$name)

# 輸出結果
print(unique_to_k)  # 輸出: 1 2
print(unique_to_j)  # 輸出: 6 7

# 找出在 df1 中但不在 df2 中的 ID
unique_to_k1 <- setdiff(k$description, j$description)

# 找出在 df2 中但不在 df1 中的 ID
unique_to_j1 <- setdiff(j$description, k$description)

# 輸出結果
print(unique_to_k1)  # 輸出: 1 2
print(unique_to_j1)  # 輸出: 6 7

# 假設有兩個數據框df1和df2
common_items <- intersect(k$description, j$description)
common_items
#改變策略，因為有些只是名字寫不同，可能是同樣的資料，所以我打算只找表格要的
pages = attributePages(ensembl_s)
pages

#格:dbSNP ID
#write.csv(listAttributes(ensembl_s), file = "ensembl_s_attr.csv")
getBM(
   attributes = "refsnp_id",
   filters = "snp_filter",
   values = "rs61749409",
   mart = ensembl_s)
#格:Gene Label
getBM(
  attributes = "associated_gene",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
#格:Gene ID
getBM(
  attributes = "ensembl_gene_name",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
#格:Chromosome
getBM(
  attributes = "chr_name",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
#格:Position
getBM(
  attributes = "chrom_start",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
#格:genome version
####ensembl_s沒有

#格:ref allele
getBM(
  attributes = "consequence_allele_string",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)

# getBM(
#   attributes = "allele",
#   filters = "snp_filter",
#   values = "rs61749409",
#   mart = ensembl_s)
# 
# getBM(
#   attributes = "allele_1",
#   filters = "snp_filter",
#   values = "rs61749409",
#   mart = ensembl_s)

#格:minor allele
#這裡反而沒資料
getBM(
  attributes = "minor_allele",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)

getBM(
  attributes = "associated_variant_risk_allele",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)

#格:MAF
#這裡反而沒資料
getBM(
  attributes = "minor_allele_freq",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
#格:Odds Ratio (OR, GWAS)
####ensembl_s沒有

#格:Population
####ensembl_s沒有
#格:Original Source (GWAS、dbSNP、Ensembl、papers)
getBM(
  attributes = "refsnp_source",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
getBM(
  attributes = "refsnp_source_description",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)
#格:p_value
#這裡反而沒資料
getBM(
  attributes = "p_value",
  filters = "snp_filter",
  values = "rs61749409",
  mart = ensembl_s)

###能取到哪些?#####
listFilters(mart = ensembl_s)

getBM(
  attributes = "refsnp_id",
  filters = "phenotype_description",
  values = "ankylosing spondylitis",
  mart = ensembl_s)
#怎就78個呢?

listFilterOptions(mart = ensembl_s, filter = "phenotype_description")
#怎麼沒有
#用來查詢就好

