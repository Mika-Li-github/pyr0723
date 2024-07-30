library(tidyverse)
getwd()
setwd("C:/Users/mikali/Desktop/gitfamily/網頁資料")
ensembl_export_Ank <- read_csv("ensembl-export_Ank.csv")
gwas_association_downloaded_EFO_0003898 <- read_tsv("gwas-association-downloaded_2024-07-30-EFO_0003898.tsv")

#先來把空的欄位去掉
summary(ensembl_export_Ank)
#name_extra, Submitter, Supporting evidence, evidence_links
ensembl_export_Ank <- ensembl_export_Ank %>% 
  select(-c('name_extra', 'Submitter', 'Supporting evidence', 'evidence_links'))

summary(gwas_association_downloaded_EFO_0003898)
#沒有

#把一定用不到的去掉
str(ensembl_export_Ank)  
#name_link, gene_links, Annotation source...11, study_links
ensembl_export_Ank <- ensembl_export_Ank %>% 
  select(-c('name_link', 'gene_links', 'Annotation source...11', 'study_links'))

str(gwas_association_downloaded_EFO_0003898)  
#DATE ADDED TO CATALOG, PUBMEDID, FIRST AUTHOR, DATE, JOURNAL, LINK, STUDY, MAPPED_TRAIT_URI
gwas_association_downloaded_EFO_0003898 <- gwas_association_downloaded_EFO_0003898 %>% 
  select(-c('DATE ADDED TO CATALOG', 'FIRST AUTHOR', 'DATE', 'JOURNAL', 'LINK', 'MAPPED_TRAIT_URI'))

#有問題SNP_ID_CURRENT(type使其跑掉)，OR or BETA(要再確認如何換算)
#去弄個大panda來
setwd("C:/Users/mikali/Desktop/gitfamily")
bigPanDa_Asso <- read_csv("bigPanDa_Asso.csv")
bigPanDa_Snps <- read_csv("bigPanDa_Snps.csv")

# 拆分 riskAlleleName
bigPanDa_Asso <- bigPanDa_Asso %>%
  mutate(riskAlleleName = as.character(riskAlleleName)) %>% # 確保是字符型
  mutate(split_pos = nchar(riskAlleleName) - 2) %>% # 計算拆分位置
  mutate(dbSNPID = substr(riskAlleleName, 1, split_pos), # 提取 dbSNPID
         riskAllele = substr(riskAlleleName, split_pos + 2, nchar(riskAlleleName))) %>% # 提取 riskAllele
  select(-split_pos) # 刪除輔助變量

# 查看結果
print(bigPanDa_Asso)

#來看看有哪些用不到的變量
str(bigPanDa_Asso)
#lastMappingDate, lastUpdateDate, riskAlleleName, loci
bigPanDa_Asso <- bigPanDa_Asso %>% 
  select(-c('lastMappingDate', 'lastUpdateDate', 'riskAlleleName','loci'))

str(bigPanDa_Snps)
#lastUpdateDate, locations, region, genomicContexts
bigPanDa_Snps <- bigPanDa_Snps %>% 
  select(-c('lastUpdateDate', 'locations', 'region', 'genomicContexts'))

#有一些變數名稱怪怪的，先看有沒有全空的值，再決定要不要刪
summary(bigPanDa_Asso)
#description...10, haplotypeSnpCount
bigPanDa_Asso <- bigPanDa_Asso %>% 
  select(-c('description...10', 'haplotypeSnpCount'))

summary(bigPanDa_Snps)
#沒有

#是不是有哪些觀測值完全重複?
# 提取所有完全重複的觀測值
duplicates <- bigPanDa_Asso %>%
  group_by(across(everything())) %>%  # 根據所有變量分組
  filter(n() > 1) %>%                  # 選擇出現次數大於1的行
  ungroup()

# 查看重複的觀測值
print(duplicates)
# bigPanDa_Asso沒有

duplicates <- bigPanDa_Snps %>%
  group_by(across(everything())) %>%  # 根據所有變量分組
  filter(n() > 1) %>%                  # 選擇出現次數大於1的行
  ungroup()

# 查看重複的觀測值
print(duplicates)
# bigPanDa_Snps有89個，我認為這代表那些SNP是一樣的

sorted_col1 <- bigPanDa_Asso %>%
  select(dbSNPID) %>%
  arrange(dbSNPID)
sorted_col1 

sorted_col2 <- bigPanDa_Snps %>%
  select(rsId) %>%
  arrange(rsId)
sorted_col2 
# 比對兩個數據框
comparison <- cbind(sorted_col1, sorted_col2)

# 找出不相等的觀測值
not_equal <- data.frame(
  sorted_col1 = sorted_col1,
  sorted_col2 = sorted_col2
) %>%
  filter(sorted_col1 != sorted_col2)

# 查看不相等的觀測值
print(not_equal)
#看來都一樣
#那要怎知是哪一筆
#依照數據，推測是不同study、不同人種的結果
#如果想知道有哪些不重複的位點，其實在excel很好做，所以我要依照study分...
#bigPanDa_Asso
#...要用AssociationId找study阿阿阿阿(python見)
duplicates <- bigPanDa_Asso %>%
  group_by(across(associationId)) %>%  # 根據associationId變量分組
  filter(n() > 1) %>%                  # 選擇出現次數大於1的行
  ungroup()
#果然是interaction他們，沒關係