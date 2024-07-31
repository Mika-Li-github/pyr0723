library(tidyverse)
getwd()
setwd("C:/Users/mikali/Desktop/gitfamily/網頁資料")
ensembl_export_Ank <- read_csv("ensembl-export_Ank.csv")
gwas_association_downloaded_EFO_0003898 <- read_tsv("gwas-association-downloaded_2024-07-30-EFO_0003898.tsv")
setwd("C:/Users/mikali/Desktop/gitfamily")
bigPanDa_Asso <- read_csv("bigPanDa_Asso.csv")
bigPanDa_SNPs <- read_csv("bigPanDa_SNPs.csv")
distinct_ensembl_export_Ank <- distinct(ensembl_export_Ank)
print(distinct_ensembl_export_Ank)
#462 to 428(沒有多算interaction?原本的表格就沒有把他們的Name(s)合在一起)

######bigPanDa_SNPs處理##########
unique_bigPanDa_SNPs <- unique(bigPanDa_SNPs)

# 拆分 riskAlleleName
bigPanDa_Asso <- bigPanDa_Asso %>%
  mutate(riskAlleleName = as.character(riskAlleleName)) %>% # 確保是字符型
  mutate(split_pos = nchar(riskAlleleName) - 2) %>% # 計算拆分位置
  mutate(dbSNPID = substr(riskAlleleName, 1, split_pos), # 提取 dbSNPID
         riskAllele = substr(riskAlleleName, split_pos + 2, nchar(riskAlleleName))) %>% # 提取 riskAllele
  select(-split_pos) # 刪除輔助變量

# 合并 A 和 B，留下共有值，并重命名列
Asso_SNPs <- bigPanDa_Asso %>%
  inner_join(unique_bigPanDa_SNPs, by = c("dbSNPID" = "rsId")) %>%
  rename("rsID" = "dbSNPID")

#來吧Gwas
# 使用separate函数分割B变量，并填充NA值
result <- gwas_association_downloaded_EFO_0003898 %>%
  separate(SNPS, into = c("dbSNP_Id", "C"), sep = " x ", fill = "right", extra = "warn", remove = FALSE)
# 查看结果
print(result)
# 互换C和B的值，前提是C不为NA
result_1 <- result %>%
  filter(!is.na(C)) %>%  # 选出C不为NA的行
  mutate(dbSNP_Id = C)  # 互换B和C的值
# 查看结果
print(result_1)
# 合并两个数据框
gwas <- bind_rows(result, result_1) %>% 
  select(-C)

#gwas的dbSNP_Id，和Asso_SNPs的rsID
# gwas_Asso_SNPs <- gwas %>%
#   inner_join(Asso_SNPs, by = c("dbSNP_Id" = "rsID")) %>%
#   rename("dbSNPId" = "rsID")
#上面不行，去看看
#都是這個HLA-B*2705，有點不想管
#但是，我們要專業
# variable_names <- names(gwas)
# print(variable_names)
# 
# variable_names_2 <- names(Asso_SNPs)
# print(variable_names_2)

# write_csv(Asso_SNPs, "Asso_SNPs.csv")
# write_csv(gwas, "gwas.csv")

#Asso_SNPs幾乎沒有作者的訊息

#先來看distinct_ensembl_export_Ank
#write_csv(distinct_ensembl_export_Ank, "distinct_ensembl_export_Ank.csv")

#distinct_ensembl_export_Ank也沒有太多作者訊息
#而且Phenotype/Disease/Trait需要正則化...
#distinct_ensembl_export_Ank$`cleaned_Phenotype/Disease/Trait` <- str_extract(distinct_ensembl_export_Ank$`Phenotype/Disease/Trait`, "(?<=\">)(.*?)(?=</a>)")
#居然搞定了，但其實不用，p_desc就是我要的
# 创建新列 PMID，提取 External reference 中 ":" 之后的内容
distinct_ensembl_export_Ank <- distinct_ensembl_export_Ank %>%
  mutate(PMID = str_extract(`External reference`, "(?<=:).*"))

write_csv(distinct_ensembl_export_Ank, "distinct_ensembl_export_Ank.csv")

#試試看用 distinct_ensembl_export_Ank 和其他欄位去對應 gwas
#distinct_ensembl_export_Ank:
#`Name(s)` = dbSNP_Id,
#`Reported gene(s)` = `REPORTED GENE(S)`, p_desc = MAPPED_TRAIT, `PMID` = PUBMEDID (改了)
#就用上面這些建立一個新欄位好了。(`Genomic location (strand)` 不等於CHR_ID和CHR_POS)

#distinct_ensembl_export_Ank
#gwas

distinct_ensembl_export_Ank$For_EGcompare <- str_c(
  str_replace_na(distinct_ensembl_export_Ank$`Name(s)`, replacement = "NA"),
  str_replace_na(distinct_ensembl_export_Ank$`Reported gene(s)`, replacement = "NA"), 
  #str_replace_na(distinct_ensembl_export_Ank$p_desc, replacement = "NA"), 
  str_replace_na(distinct_ensembl_export_Ank$`PMID`, replacement = "NA"), 
  sep = " ")

gwas$For_EGcompare <- str_c(
  str_replace_na(gwas$dbSNP_Id, replacement = "NA"),
  str_replace_na(gwas$`REPORTED GENE(S)`, replacement = "NA"), 
  #str_replace_na(gwas$MAPPED_TRAIT, replacement = "NA"), 
  str_replace_na(gwas$PUBMEDID, replacement = "NA"), 
  sep = " ")

result <- inner_join(distinct_ensembl_export_Ank, gwas, by = "For_EGcompare")
print(result)













