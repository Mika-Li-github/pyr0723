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
len(snps.variants.at[0,'genomicContexts'])
