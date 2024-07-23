# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 15:53:25 2024

@author: mikali
"""
import pandas as pd


from pandasgwas.get_variants import get_variants_by_study_id
variants = get_variants_by_study_id('GCST002318')

'''
from pandasgwas.get_studies import get_studies_by_study_id
studies = get_studies_by_study_id('GCST002318')
'''

from pandasgwas.get_associations import get_associations_by_efo_trait
associations = get_associations_by_efo_trait('ankylosing spondylitis')

#需要:associations、variants
type(associations.raw_data)
associations.raw_data[0]['orPerCopyNum']
associations.raw_data[0]['pvalue']

data = {
    'orPerCopyNum': associations.raw_data[0]['orPerCopyNum'],
    'pvalue': associations.raw_data[0]['pvalue']
}

df = pd.DataFrame([data])
print(df)


associations.loci

indices = [index for index, item in enumerate(associations.raw_data) if len(item['loci']) > 1]

print(indices)
#[349, 350]是SNP x SNP interaction，算兩個?，associationse共326個SNP
#應該是看strongest_risk_alleles，用associationid挖資料
#先把rsID找出來
associations.raw_data[0]["loci"]
associations.raw_data[349]["loci"][0]
#果然是看strongest_risk_alleles，再用associationid挖資料

# 假设 associations.strongest_risk_alleles 是一个 DataFrame
# 找出 associationId 重复的索引
duplicate_indices = associations.strongest_risk_alleles[associations.strongest_risk_alleles.duplicated(subset='associationId', keep=False)].index

print(duplicate_indices)
#重複的有Index([349, 350, 351, 352], dtype='int64')，可行

# 假设 associations.strongest_risk_alleles 是一个 DataFrame
# 提取 associationId 和 riskAlleleName 列
new_df = associations.strongest_risk_alleles[['associationId', 'riskAlleleName']].copy()
# 输出新的 DataFrame
print(new_df)

# 从 associations 中选择所需列
associations_subset = associations.associations[['associationId', 'orPerCopyNum', 'pvalue']]

# 合并 DataFrame
merged_df = pd.merge(new_df, associations_subset, on='associationId', how='left')

# 输出合并后的 DataFrame
print(merged_df)

associations.associations.at[425,"orPerCopyNum"]
