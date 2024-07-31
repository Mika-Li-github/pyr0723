library(tidyverse)
getwd()
setwd("C:/Users/mikali/Desktop/gitfamily/網頁資料")
ensembl_export_Ank <- read_csv("ensembl-export_Ank.csv")
gwas_association_downloaded_EFO_0003898 <- read_tsv("gwas-association-downloaded_2024-07-30-EFO_0003898.tsv")
setwd("C:/Users/mikali/Desktop/gitfamily")
bigPanDa_Asso <- read_csv("bigPanDa_Asso.csv")
bigPanDa_SNPs <- read_csv("bigPanDa_SNPs.csv")

#來檢查數據
head(ensembl_export_Ank)
str(ensembl_export_Ank)
summary(ensembl_export_Ank)
#格 dbSNP ID: Name(s)

head(gwas_association_downloaded_EFO_0003898)
str(gwas_association_downloaded_EFO_0003898)
summary(gwas_association_downloaded_EFO_0003898)
#格 dbSNP ID: SNPS
#376 377有問題，應該先分開

head(bigPanDa_Asso)
str(bigPanDa_Asso)
summary(bigPanDa_Asso)
#格 dbSNP ID: riskAlleleName
# 拆分 riskAlleleName
bigPanDa_Asso <- bigPanDa_Asso %>%
  mutate(riskAlleleName = as.character(riskAlleleName)) %>% # 確保是字符型
  mutate(split_pos = nchar(riskAlleleName) - 2) %>% # 計算拆分位置
  mutate(dbSNPID = substr(riskAlleleName, 1, split_pos), # 提取 dbSNPID
         riskAllele = substr(riskAlleleName, split_pos + 2, nchar(riskAlleleName))) %>% # 提取 riskAllele
  select(-split_pos) # 刪除輔助變量

head(bigPanDa_SNPs)
str(bigPanDa_SNPs)
summary(bigPanDa_SNPs)
#格 dbSNP ID: rsId


# 使用separate函数分割B变量，并填充NA值
result <- gwas_association_downloaded_EFO_0003898 %>%
  separate(SNPS, into = c("dbSNP_Id", "C"), sep = " x ", fill = "right", extra = "drop", remove = FALSE)

# 查看结果
print(result)


# 互换C和B的值，前提是C不为NA
result_1 <- result %>%
  filter(!is.na(C)) %>%  # 选出C不为NA的行
  mutate(dbSNP_Id = C)  # 互换B和C的值

# 查看结果
print(result_1)

# 合并两个数据框
gwas_association_downloaded_EFO_0003898_combined_result <- bind_rows(result, result_1) %>% 
  select(-C)

#看來，是時候合併數據了
######先來看看重複的值有哪些異同######

# 查找列A中重复的观测值及其完整信息
result <- ensembl_export_Ank %>%
  group_by(`Name(s)`) %>%                # 按A列分组
  filter(n() > 1) %>%            # 筛选出重复的观测值
  ungroup()                      # 取消分组

# 查看结果
print(result)

# 查找重复的整列观测值
duplicates <- ensembl_export_Ank[duplicated(ensembl_export_Ank) | duplicated(ensembl_export_Ank, fromLast = TRUE), ]

# 查看结果
print(duplicates)
#移除吧
# 找到重复的整列观测值并保留唯一记录
duplicates <- ensembl_export_Ank[duplicated(ensembl_export_Ank) | duplicated(ensembl_export_Ank, fromLast = TRUE), ] %>%
  distinct()

# 查看结果
print(duplicates)
#33個
distinct_data <- distinct(ensembl_export_Ank)
print(distinct_data)
#462 to 428(沒有多算interaction?原本的表格就沒有把他們的Name(s)合在一起)


# 找到重复的整列观测值并保留唯一记录
duplicates <- gwas_association_downloaded_EFO_0003898_combined_result[duplicated(gwas_association_downloaded_EFO_0003898_combined_result) | duplicated(gwas_association_downloaded_EFO_0003898_combined_result, fromLast = TRUE), ] %>%
  distinct()

# 查看结果
print(duplicates)
#0個

# 找到重复的整列观测值并保留唯一记录
duplicates <- bigPanDa_Asso[duplicated(bigPanDa_Asso) | duplicated(bigPanDa_Asso, fromLast = TRUE), ] %>%
  distinct()

# 查看结果
print(duplicates)
#0個

# 找到重复的整列观测值并保留唯一记录
duplicates <- bigPanDa_SNPs[duplicated(bigPanDa_SNPs) | duplicated(bigPanDa_SNPs, fromLast = TRUE), ] %>%
  distinct()

# 查看结果
print(duplicates)
#38個

distinct_data <- distinct(bigPanDa_SNPs)
print(distinct_data)
#426 to 375 (算了interaction)

#來看看有沒有只出現在bigPanDa_SNPs的rsId
# 使用 anti_join 找出只在 A 表格中的 B 变量
only_in_bigPanDa_SNPs <- anti_join(bigPanDa_SNPs, bigPanDa_Asso, by = c("rsId" = "dbSNPID"))
print(only_in_bigPanDa_SNPs)
#沒有，那反過來呢?

only_in_bigPanDa_Asso <- anti_join(bigPanDa_Asso, bigPanDa_SNPs, by = c("dbSNPID" = "rsId"))
print(only_in_bigPanDa_Asso)
#也沒有
common_values <- inner_join(bigPanDa_Asso, bigPanDa_SNPs, by = c("dbSNPID" = "rsId"))
print(common_values)
unique_data <- unique(common_values)
#也就是兩者的rsId完全重複

#來比較bigPanDa_SNPs和ensembl_export_Ank的差別
only_in_bigPanDa_SNPs_vs_E <- anti_join(bigPanDa_SNPs, ensembl_export_Ank, by = c("rsId" = "Name(s)"))
print(only_in_bigPanDa_SNPs_vs_E)
#23個

#來比較bigPanDa_SNPs和gwas_association_downloaded_EFO_0003898_combined_result的差別
only_in_bigPanDa_SNPs_vs_G <- anti_join(bigPanDa_SNPs, gwas_association_downloaded_EFO_0003898_combined_result, by = c("rsId" = "dbSNP_Id"))
print(only_in_bigPanDa_SNPs_vs_G)
#0個，不意外

common_values_G <- inner_join(bigPanDa_SNPs, gwas_association_downloaded_EFO_0003898_combined_result, by = c("rsId" = "dbSNP_Id"))
print(common_values_G)
unique_data <- unique(common_values_G)
print(unique_data)
#有一樣的，也就是兩者完全重複，畢竟都是GWAS來的。
#先把三個一樣的合在一起
