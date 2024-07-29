library(tidyverse)
library(biomaRt)
library(BiocManager)

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
filters <- listFilters(ensembl)
print(filters)
listFilterOptions(mart = ensembl, filter = "phenotype_description")
#沒有Ankylosing spondylitis，也沒有其他疾病，卻可用。
listAttributes(mart = ensembl, what = "name", page = "snp")
as_bio_attr <- listAttributes(mart = ensembl, page = "snp")
write.csv(as_bio_attr, file = "as_bio_attr")
as_attributes <- c(
  "refsnp_id",
  "refsnp_source",
  "chr_name",
  "chrom_start",
  "chrom_end",
  "allele",
  "allele_1",
  "minor_allele",
  "minor_allele_freq",
  "clinical_significance", #怎麼變103個
  "associated_gene", #怎麼變104個
  #"phenotype_description",
  "associated_variant_risk_allele",
  "p_value"#,變137個，有很多重複的snp項?
  #"set_name"#變兩千多個，有很多重複的項?
  
)

results <- getBM(attributes = as_attributes,
                 filters = 'phenotype_description',
                 values = 'Ankylosing spondylitis',
                 mart = ensembl)
print(results)

#相關基因太少，重複項的問題，我打算只用來做確認

#######來整理其他檔案########
#先看有哪些snp(唯一rsID)
#叫SNPS，在gwas-association-downloaded_2024-07-26-EFO_0003898
#叫Name(s)，在ensembl-export_Loci associated with Ankylosing spondylitis
#叫dbSNPID，在asPandasGWAS
#我先去整理那些檔案

#asPandasGWAS
#Variant: imm_16_28525386，Variant: HLA-B*2705, Variant: chr18:14723700,Variant: rs67025039
#會刪除上面這些Variant does not map to the genome

#ensembl-export_Loci associated with Ankylosing spondylitis
#清掉了chromosome標了可能是載體的

#gwas-association-downloaded_2024-07-26-EFO_0003898
#chr18:14723700-? HLA-B*2705-? rs67025039-? imm_16_28525386-A
##注意，有這種的我不刪
#5 x 6	96788627 x 31377139
#5 x 6	96812030 x 31377139

#在工作表「不想重複」，我把snp交互作用的分開算成兩個位點
#當前有382個rsId

######382SNPs#######
#先把那些SNP召喚進來吧
#getwd()
#setwd("C:/Users/mikali/Desktop/gitfamily")
asPandasGWAS <- read_csv("asPandasGWAS.csv")
ensembl_export_LociAssociatedWithAnkylosingSpondylitis <- read_csv("ensembl-export_Loci associated with Ankylosing spondylitis.csv")
gwas_association_downloaded_2024_07_26_EFO_0003898 <- read_csv("gwas-association-downloaded_2024-07-26-EFO_0003898.csv")


# 找出重複的 dbSNPID
duplicates <- asPandasGWAS %>%
  group_by(dbSNPID) %>%
  summarise(count = n()) %>%
  filter(count > 1)

# 提取所有重複的觀測值
result <- asPandasGWAS %>%
  filter(dbSNPID %in% duplicates$dbSNPID)

# 查看結果
print(result, n = 85)

# 按照 orPerCopyNum 進行降序排列
sorted_data <- result %>%
  arrange(desc(dbSNPID))

# 查看結果
print(sorted_data, n = 85)



# 找出重複的 Name(s)
duplicates <- ensembl_export_LociAssociatedWithAnkylosingSpondylitis %>%
  group_by(`Name(s)`) %>%
  summarise(count = n()) %>%
  filter(count > 1)

# 提取所有重複的觀測值
result <- ensembl_export_LociAssociatedWithAnkylosingSpondylitis %>%
  filter(`Name(s)` %in% duplicates$`Name(s)`)

# 查看結果
print(result, n = 62)

# 按照 orPerCopyNum 進行降序排列
sorted_data <- result %>%
  arrange(desc(`Name(s)`))

# 查看結果
print(sorted_data, n = 62)




# 找出重複的 SNPS
duplicates <- gwas_association_downloaded_2024_07_26_EFO_0003898 %>%
  group_by(SNPS) %>%
  summarise(count = n()) %>%
  filter(count > 1)

# 提取所有重複的觀測值
result <- gwas_association_downloaded_2024_07_26_EFO_0003898 %>%
  filter(SNPS %in% duplicates$SNPS)

# 查看結果
print(result,SNPS)

# 按照 orPerCopyNum 進行降序排列
sorted_data <- result %>%
  arrange(desc(SNPS))

# 查看結果
print(sorted_data, n = 81)
#我覺得pandas可以再補完整一點