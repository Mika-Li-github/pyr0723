# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:02:07 2024

@author: mikali
"""
import os

curr_path = os.getcwd()
print(curr_path)
os.chdir('C:/Users/mikali/Desktop/gitfamily/bigPanda')

import pandas as pd

bigPanDa_Asso = pd.read_csv('bigPanDa_Asso.csv')
print(bigPanDa_Asso)

bigPanDa_Snps = pd.read_csv('bigPanDa_Snps.csv')
print(bigPanDa_Snps)

os.chdir('C:/Users/mikali/Desktop/gitfamily/網頁資料/0801Vectordeleted')

ensembl_export_Ank = pd.read_csv('ensembl-export_Ank.csv')
print(ensembl_export_Ank)

'''
import csv
import chardet

# 自動檢測 TSV 檔案的編碼格式
with open('gwas-association-downloaded_2024-07-30-EFO_0003898.tsv', 'rb') as tsvfile:
    result = chardet.detect(tsvfile.read())
    encoding = result['encoding']

# 使用檢測到的編碼格式讀取 TSV 檔案    
with open('gwas-association-downloaded_2024-07-30-EFO_0003898.tsv', 'r', encoding=encoding) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    rows = list(reader)
# 寫入 CSV 檔案
with open('output.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(rows)
'''

gwas = pd.read_csv('output.csv')
print(gwas)

# 提取三個表格中的該變數並合併
combined_values = pd.concat([
    bigPanDa_Snps['rsId'],
    ensembl_export_Ank['Name(s)'],
    gwas['SNPS']
])

# 選出唯一值
unique_values = combined_values.unique()

# 輸出唯一值
print(f"唯一值: {unique_values}")
len(unique_values)

# 使用 str.rsplit() 方法從最右邊的 '-' 切開變數
bigPanDa_Asso[['rsId', 'riskAllele']] = bigPanDa_Asso['riskAlleleName'].str.rsplit('-', n=1, expand=True)

# 輸出結果
print(bigPanDa_Asso)


ensembl_export_Ank.rename(columns={'Name(s)': 'rsId'}, inplace=True)

gwas.rename(columns={'SNPS': 'rsId'}, inplace=True)

rsids = unique_values

import pandas as pd


# 事先決定好要查詢的 rsId 陣列
query_rsIds = unique_values

import pandas as pd

def filter_first_occurrences(df, query_rsIds):
    """
    根據給定的 rsId 陣列，過濾出資料表中每個 rsId 第一次出現的觀測值。

    參數:
    df (pd.DataFrame): 要查詢的資料表。
    query_rsIds (list): 要查詢的 rsId 陣列。

    返回:
    pd.DataFrame: 包含每個 rsId 第一次出現的觀測值的新表格。
    """
    # 過濾出在查詢陣列中的 rsId
    filtered_df = df[df['rsId'].isin(query_rsIds)]

    # 找到每個 rsId 的第一次出現
    new_table = filtered_df.drop_duplicates(subset='rsId', keep='first')

    return new_table

# 使用函式獲取新的表格
new_table_bigPanDa_Asso = filter_first_occurrences(bigPanDa_Asso, query_rsIds)

# 顯示新表格
print(new_table_bigPanDa_Asso)

new_table_bigPanDa_Snps = filter_first_occurrences(bigPanDa_Snps, query_rsIds)

new_table_ensembl_export_Ank = filter_first_occurrences(ensembl_export_Ank, query_rsIds)

new_table_gwas = filter_first_occurrences(gwas, query_rsIds)


import pandas as pd

# 使用 merge 合併表格
merged_data = pd.merge(new_table_bigPanDa_Asso, new_table_bigPanDa_Snps, on='rsId', how='outer')
merged_data = pd.merge(merged_data, new_table_ensembl_export_Ank, on='rsId', how='outer')
merged_data = pd.merge(merged_data, new_table_gwas, on='rsId', how='outer')

# 顯示合併後的表格
print(merged_data)

#搞...搞定了?
#我想選出
new_df = merged_data[['rsId', 'CHR_ID', 'CHR_POS', 'chromosomeName', 'chromosomePosition', 'Genomic location (strand)']]

#想檢查這些變數的意義


# 選擇除了 'rsId' 列外，其他列完全沒有 NA 值的觀測值
new_df_noNA = new_df.dropna(subset=new_df.columns.difference(['rsId']))

# 顯示新的 DataFrame
print(new_df_noNA)


# 提取數字，如果是 'X' 則保留原值，並將所有數字轉換為字符串
new_df_noNA['提取的數字'] = new_df_noNA['chromosomeName'].apply(lambda x: str(int(eval(x)[0])) if x != "['X']" else 'X')

# 顯示新的 DataFrame
print(new_df_noNA)


# 比較 '提取的數字' 和 'CHR_ID' 欄位是否完全相等
are_equal = new_df_noNA['提取的數字'].equals(new_df_noNA['CHR_ID'])

# 顯示比較結果
print(f"兩個欄位是否完全相等: {are_equal}")
#兩個欄位是否完全相等: True
#可以輸出表格了

new_df.to_csv('PandaAndWebSNP.csv', index=False, encoding='utf-8-sig')
