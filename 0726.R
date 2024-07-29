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
#用來查詢就好?

## Use the Ensembl human genes dataset

ensembl_s1 <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")

## we can search for the name of a filter we're interested in e.g. 'phenotype'
## we need to use the name of the filter in the next function
searchFilters(ensembl_s1, pattern = "phenotype")

## list all the options available to the 'phenotype_source' filter
listFilterOptions(mart = ensembl_s1, filter = "phenotype_source")

## search the 'phenotype_description' filter for the term 'as'
searchFilterOptions(mart = ensembl_s1,
                    filter = "phenotype_description",
                    pattern = "as")
#還是不行
listFilters(ensembl_s)

listFilters(ensembl_g)
searchFilters(ensembl_g, pattern = "pheno")
listFilterOptions(mart = ensembl_g, filter = "phenotype_description")
searchFilterOptions(mart = ensembl_g, filter = "phenotype_description", pattern = "ank")
#有耶
searchFilterOptions(mart = ensembl_s, filter = "phenotype_description", pattern = "ank")

listAttributes(ensembl_s)

library(biomaRt)

# 连接到Ensembl数据库
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# 检查可用的过滤器
filters <- listFilters(ensembl)
print(filters)

# 列出phenotype_description过滤器的选项
options <- listFilterOptions(ensembl, "phenotype_description")
print(options)

# 使用phenotype_description过滤器查询
results <- getBM(attributes = c('refsnp_id', 'phenotype_description', 'phenotype_source'),
                 filters = 'phenotype_description',
                 values = 'Ankylosing spondylitis',
                 mart = ensembl)

print(results)
#還是那78個，好窩
listAttributes(mart = ensembl)
listAttributes(mart = ensembl_s1)

searchAttributes(mart = ensembl, pattern = "phe")
searchAttributes(mart = ensembl_s1, pattern = "phe")

results <- getBM(attributes = c('refsnp_id', 'phenotype_description', 'phenotype_name'),
                 filters = 'phenotype_description',
                 values = 'Ankylosing spondylitis',
                 mart = ensembl)

print(results)

results <- getBM(attributes = c('refsnp_id', 'phenotype_description', 'phenotype_name','validated'),
                 filters = 'phenotype_description',
                 values = 'Ankylosing spondylitis',
                 mart = ensembl_s1)

print(results)
#雖然都沒有寫來源，我還是想看是哪78個
resultsS <- getBM(attributes = c('refsnp_id', 'phenotype_description', 'phenotype_name'),
                 filters = 'phenotype_description',
                 values = 'Ankylosing spondylitis',
                 mart = ensembl_s1)

print(resultsS)
#write.csv(resultsS, file="as78biomaRt.csv")
#除了rs55688811，其他都和pandasGWAS結果重複

#先把想要的變數找出來，再看要怎麼比對