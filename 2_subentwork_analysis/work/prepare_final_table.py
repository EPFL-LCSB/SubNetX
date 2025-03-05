#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 16:22:43 2024

@author: omid
"""

import os
import numpy as np
import pandas as pd


from subnetx.io.json import load_json_model, save_json_model


from eval_netws import target_list
from enumerate_subnets import path_load, path_mod, path_enum_pthws
# from find_minimal_subnets import product_threshold
from extract_min_pathways import path_save
from ext_ids_min_path import new_path_save

# User input
MAX_LENGTH = 10
final_path_save = '{}/final_table.csv'
THRESHOLDS = [0.1, 0.5, 1]


if __name__ == "__main__":
    
    for target in target_list:
        ### preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        reactions = pd.read_csv(the_path + '/reactions.tsv', sep='\t', index_col = 0)
        compounds = pd.read_csv(the_path + '/compounds.tsv', sep='\t', index_col = 0)
        
        blocks = []
        for product_threshold in THRESHOLDS:
            # create an empty block
            this_block = pd.DataFrame()
            
            data_1 = pd.read_csv(path_enum_pthws.format(target), index_col = 0)
            data_2 = pd.read_csv(path_save.format(the_path, product_threshold), index_col = 0)
            data_3 = pd.read_csv(new_path_save.format(the_path, product_threshold), index_col = 0)
            
            # the maximum flux yield to normalize other yields
            max_yield = pd.concat([data_1]*len(data_2), ignore_index=True)
            this_block['fraction maximum yield'] = \
                data_2['rel_prod_flux']/(max_yield['product_flux'] / max_yield['available_substrate'])
            this_block['number of heterologous reactions'] = data_2['hetero_length']
            
            # extracting information for reactions
            rxn_block = pd.DataFrame(index=range(len(data_2)),columns=range(MAX_LENGTH))
            for ind, row in data_3.iterrows():
                rxn_ids = row['reactions'].split(' // ')
                for index, id_ in enumerate(rxn_ids):
                     rxn = reactions.loc[int(id_)]['R_PR_STOICH']
                     # replace the LCSB IDs with human understandable names
                     helper  = rxn.split(' ')
                     for i, x in enumerate(helper):
                         try:
                             y = compounds.loc[int(x)]['M_PR_NAME']
                             helper[i] = y
                         except KeyError:
                             pass
                         except ValueError:
                             pass
                     rxn_block[index][ind] = ' '.join(helper)
            this_block_thres = pd.concat([this_block, rxn_block], axis=1)
            blocks += [this_block_thres]
        final_table = pd.concat(blocks, axis = 0)
        final_table = final_table.drop_duplicates()
        final_table.to_csv(final_path_save.format(the_path))                     
            
