# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 11:42:11 2024

@author: mikali
"""

from pandasgwas import get_child_efo
child_list1=get_child_efo('EFO_0003898')
#MONDO_0020655 是啥

from pandasgwas.get_traits import get_traits_by_efo_id
traits = get_traits_by_efo_id(efo_id='MONDO_0020655')
#沒有traits?那我可不管了

from pandasgwas.get_traits import get_traits_by_efo_id
traits = get_traits_by_efo_id(efo_id='EFO_0003898')
from pandasgwas.get_traits import get_traits_by_efo_trait
traits = get_traits_by_efo_trait('ankylosing spondylitis')
#EFO_0003898和ankylosing spondylitis 互相唯一對應

#這次嘗試把所有能取得的資訊都先納入表格
from pandasgwas.get_associations import get_associations_by_efo_id
associations = get_associations_by_efo_id('EFO_0003898')

from pandasgwas.get_studies import get_studies_by_efo_id
studies = get_studies_by_efo_id('EFO_0003898')

from pandasgwas.get_variants import get_variants_by_efo_id
snps = get_variants_by_efo_id('EFO_0003898')
#以這個表為主，加入其他資訊

import numpy as np

# 定義一個函數來提取 key
def extract_keys(locations, key_name):
    if locations is None:
        return np.nan  # 如果 locations 為 None，返回 NA
    else:
        return [d.get(key_name, np.nan) for d in locations]  # 提取 key

# 提取 key1 和 key2
bigPanda = snps.variants
bigPanda['chromosomeName'] = bigPanda['locations'].apply(lambda loc: extract_keys(loc, 'chromosomeName'))
bigPanda['chromosomePosition'] = bigPanda['locations'].apply(lambda loc: extract_keys(loc, 'chromosomePosition'))
bigPanda['region'] = bigPanda['locations'].apply(lambda loc: extract_keys(loc, 'region'))
bigPanda['region_name'] = bigPanda['region'].apply(lambda loc: extract_keys(loc, 'name'))


# 定義一個函數來過濾符合條件的 dict
def filter_genomic_contexts(contexts):
    if contexts is None:
        return np.nan  # 如果 contexts 為 None，返回 NA
    else:
        # 過濾符合條件的 dict
        filtered_dicts = [d for d in contexts if d.get('source') == 'Ensembl' and d.get('distance') == 0]
        return filtered_dicts if filtered_dicts else np.nan  # 如果沒有符合的 dict，返回 NA

# 應用過濾函數
bigPanda['filtered_genomicContexts'] = bigPanda['genomicContexts'].apply(filter_genomic_contexts)

# 顯示結果
print(bigPanda[['genomicContexts', 'filtered_genomicContexts']])



# 定義一個函數來計算 list 中 dict 的個數
def count_dicts(lst):
    if isinstance(lst, list):
        return len(lst)
    else:
        return 0  # 如果不是 list，返回 0

# 計算每個元素中 dict 的個數
bigPanda['num_genes'] = bigPanda['filtered_genomicContexts'].apply(count_dicts)

# 顯示結果
print(bigPanda)

#snps大致上搞定了，現在要來看看association
bigPandA = associations.associations

# 使用 explode 方法將 loci 列中的 list 展開為多行
bigPandA = bigPandA.explode('loci')

import pandas as pd

# 假設 bigPandA 是您的 DataFrame
# 假設 loci 是包含字典的變量

# 創建一個新的 DataFrame 用於存儲提取的鍵值對
new_df = pd.DataFrame()

# 遍歷每一行
for index, row in bigPandA.iterrows():
    # 獲取每一行的 loci 值（字典）
    loci_dict = row['loci']
    
    # 確保 loci_dict 是字典
    if isinstance(loci_dict, dict):
        # 將字典轉換為 DataFrame，鍵作為列名，值作為觀測值
        temp_df = pd.DataFrame([loci_dict])  # 將字典放入列表中以創建 DataFrame
        new_df = pd.concat([new_df, temp_df], ignore_index=True)

# 將新的 DataFrame 添加到原始 DataFrame 中
bigPandA = pd.concat([bigPandA.reset_index(drop=True), new_df.reset_index(drop=True)], axis=1)

# 顯示結果
print(bigPandA)

import pandas as pd

# 假設 bigPandA 是您的 DataFrame
# 假設 strongestRiskAlleles 是包含字典列表的變量

# 創建一個新的 DataFrame 用於存儲提取的鍵值對
new_df = pd.DataFrame()

# 遍歷每一行
for index, row in bigPandA.iterrows():
    # 獲取每一行的 strongestRiskAlleles 值（列表）
    alleles_list = row['strongestRiskAlleles']
    
    # 確保 alleles_list 是列表
    if isinstance(alleles_list, list):
        # 遍歷列表中的每個字典
        for allele_dict in alleles_list:
            if isinstance(allele_dict, dict):
                # 將字典的鍵值對展平到新的 DataFrame
                temp_df = pd.DataFrame([allele_dict])  # 將字典放入列表中以創建 DataFrame
                new_df = pd.concat([new_df, temp_df], ignore_index=True)

# 將新的 DataFrame 添加到原始 DataFrame 中
bigPandA = pd.concat([bigPandA.reset_index(drop=True), new_df.reset_index(drop=True)], axis=1)

# 顯示結果
print(bigPandA)


#先來去除一些不需要的變數吧
#去R操作好了
bigPanDa_Asso = bigPandA
bigPanDa_Snps = bigPanda


import os

# 獲取當前工作目錄
current_directory = os.getcwd()

# 打印當前工作目錄
print("Current Directory:", current_directory)

os.chdir('C:/Users/mikali/Desktop/gitfamily')

files = os.listdir(current_directory)
print("Files in Current Directory:", files)

import pandas as pd
# 假設 bigPandA 是您要寫入 CSV 的 DataFrame
# 將 DataFrame 寫入 CSV 檔案
bigPanDa_Asso.to_csv('bigPanDa_Asso.csv', index=False, encoding='utf-8')
bigPanDa_Snps.to_csv('bigPanDa_Snps.csv', index=False, encoding='utf-8')

