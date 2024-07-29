# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:01:03 2024

@author: mikali
"""
#pip install pandasgwas

#Trait: EFO_0003898
from pandasgwas.get_associations import get_associations_by_efo_id
associations = get_associations_by_efo_id('EFO_0003898')

#格:dbSNP ID
associations.strongest_risk_alleles.at[0, 'riskAlleleName']
#另一個選擇
from pandasgwas.get_variants import get_variants_by_efo_id
snps = get_variants_by_efo_id('EFO_0003898')

snps.variants.at[0,'rsId']

#格:Gene Label(associations有基因沒距離，snps有，建議一起取出)
#distance=0, source=Ensemble
#那我要改變一下，在做asSnp的時候，gene不要把key提取出來變欄位，就直接一欄叫gene
import pandas as pd

# 初始化一个空列表来存储每个索引的字典
data = []

# 遍历 snps.raw_data
for entry in snps.raw_data:
    # 提取 rsId
    rsId = entry['rsId']
    
    # 检查 locations 列表是否为空
    if entry['locations']:
        # 获取第一个 location dict
        location = entry['locations'][0]
        
        # 提取 chromosomeName 和 chromosomePosition
        chromosomeName = location['chromosomeName']
        chromosomePosition = location['chromosomePosition']
    else:
        # 如果 locations 为空,设置默认值
        chromosomeName = None
        chromosomePosition = None
    
    # 提取所有符合条件的 gene
    genes = []  # 用于存储符合条件的 gene
    for context in entry['genomicContexts']:
        if context['distance'] == 0 and context['source'] == 'Ensembl':
            if 'gene' in context:
                # 只提取所需的字段
                gene_info = {
                    'geneName': context['gene'].get('geneName'),
                    'entrezGeneIds': context['gene'].get('entrezGeneIds'),
                    'ensemblGeneIds': context['gene'].get('ensemblGeneIds')
                }
                genes.append(gene_info)  # 将 gene 信息字典添加到列表中
    
    # 创建记录字典
    record = {
        'rsId': rsId,
        'chromosomeName': chromosomeName,
        'chromosomePosition': chromosomePosition,
        'gene': genes
    }
    
    # 将当前记录添加到数据列表
    data.append(record)

# 创建 DataFrame
asSnp = pd.DataFrame(data)

# 显示结果
print(asSnp)



#好了，snp搞定，現在來看association能提供甚麼資料

import pandas as pd

# 初始化一个空列表来存储每个索引的字典
data = []

# 遍历 associations.raw_data
for entry in associations.raw_data:
    # 提取所需的字段
    pvalue = entry['pvalue']
    orPerCopyNum = entry['orPerCopyNum']
    associationId = entry['associationId']
    
    # 检查 loci 是否为空
    if entry['loci']:  # 如果 loci 列表不为空
        for locus in entry['loci']:
            strongestRiskAlleles = locus.get('strongestRiskAlleles', [])
            if strongestRiskAlleles:  # 检查 strongestRiskAlleles 是否为空
                for allele in strongestRiskAlleles:
                    riskAlleleName = allele.get('riskAlleleName', None)  # 获取 riskAlleleName
                    
                    # 创建记录字典
                    record = {
                        'pvalue': pvalue,
                        'orPerCopyNum': orPerCopyNum,
                        'riskAlleleName': riskAlleleName,
                        'associationId': associationId
                    }
                    
                    # 将当前记录添加到数据列表
                    data.append(record)

# 创建 DataFrame
asSoc = pd.DataFrame(data)

# 显示结果
print(asSoc)

# 查找重复的 associationId
duplicate_association_ids = asSoc[asSoc.duplicated(subset='associationId', keep=False)]

# 获取重复的索引
duplicate_indices = duplicate_association_ids.index.tolist()

# 显示结果
print("重复的 associationId 的索引:", duplicate_indices)

# 使用字符串操作从右侧分割 riskAlleleName
split_risk_alleles = asSoc['riskAlleleName'].str.rsplit('-', n=1, expand=True)

# 将分割结果赋值给新的列
asSoc['riskRsId'] = split_risk_alleles[0]
asSoc['riskAllele'] = split_risk_alleles[1]

# 删除原来的 riskAlleleName 列
asSoc = asSoc.drop(columns=['riskAlleleName'])

# 显示更新后的 DataFrame
print(asSoc)

#penpinappleapplepen
# 检查 asSoc 的 riskRsId 是否能完全与 asSnp 的 rsId 对应
if all(asSoc['riskRsId'].isin(asSnp['rsId'])):
    # 合并 DataFrame
    asPandasGWAS = pd.merge(asSoc, asSnp, left_on='riskRsId', right_on='rsId', how='inner')
    
    # 重命名 rsId 列为 dbSNPID
    asPandasGWAS.rename(columns={'rsId': 'dbSNPID'}, inplace=True)
    
    # 删除 riskRsId 列
    asPandasGWAS.drop(columns=['riskRsId'], inplace=True)
    # 假设 asPandasGWAS 是您的 DataFrame
    # 将 'gene' 列转换为字符串
    asPandasGWAS['gene'] = asPandasGWAS['gene'].apply(lambda x: str(x) if isinstance(x, list) else x)

    # 确保合并后的 DataFrame 不包含重复行
    asPandasGWAS = asPandasGWAS.drop_duplicates()
    
    # 显示合并后的 DataFrame
    print(asPandasGWAS)
else:
    print("asSoc 的 riskRsId 并不完全能与 asSnp 的 rsId 对应。")
    

    
# 将 asPandasGWAS 保存为 CSV 文件
asPandasGWAS.to_csv('asPandasGWAS.csv', index=False)

# 初始化一个空列表来存储每个索引的基因数量
gene_counts = []

# 遍历 snps.raw_data
for entry in snps.raw_data:
    gene_count = 0
    
    # 检查 genomicContexts
    for context in entry['genomicContexts']:
        if context['distance'] == 0 and context['source'] == 'Ensembl':
            if 'gene' in context:
                gene_count += 1
    
    # 将当前索引的基因数量添加到列表中
    gene_counts.append(gene_count)

# 检查是否所有索引都只有一个基因
if all(count == 1 for count in gene_counts):
    print("每个索引都只有一个符合条件的基因。")
else:
    print("有些索引有多个符合条件的基因。")
    print("以下是每个索引的基因数量:", gene_counts)

