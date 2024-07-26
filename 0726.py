# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:01:03 2024

@author: mikali
"""
#Trait: EFO_0003898
from pandasgwas.get_associations import get_associations_by_efo_id
associations = get_associations_by_efo_id('EFO_0003898')

from pandasgwas.get_associations import get_associations_by_efo_trait
associations = get_associations_by_efo_trait('Ankylosing spondylitis')
