#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 16:48:46 2021

@author: omid
"""
import os
import pandas as pd
import numpy as np

from tqdm import tqdm

from enumerate_subnets import target_list, path_mod
from find_minimal_subnets import path_enum_pthws, product_threshold

BIN_VAR_PREFIX = 'BFUSE_'
path_save = '{}/pthws_info_{}_min.csv'

if __name__ == "__main__":
    
    for target in target_list:
        ### preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        data = pd.read_csv(path_enum_pthws.format(target, product_threshold))
        pathways = dict()
        for ind, row in tqdm(data.iterrows(), desc= 'scanning the pathways'):
            pthw_dict = {'rel_prod_flux' : \
                         row['product_flux']/row['available_substrate']}
            pthw_rxns = []
            for id_, bin_var in row.iteritems():
                if BIN_VAR_PREFIX in id_ and np.isclose(bin_var,1):
                    pthw_rxns.append(id_.replace(BIN_VAR_PREFIX,''))
            
            pthw_dict['hetero_length'] = len(pthw_rxns)
            pthw_dict['reactions'] = pthw_rxns
            pathways[ind] = pthw_dict
        pd.DataFrame.from_dict(pathways,
                               orient='index').to_csv(
                                   path_save.format(the_path,
                                                    product_threshold))
        