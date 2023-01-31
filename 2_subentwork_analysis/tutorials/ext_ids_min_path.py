#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:52:52 2022

@author: omid
"""


import os
import numpy as np
import pandas as pd

from tqdm import tqdm

from pytfa.optim.variables import ForwardBackwardUseVariable
from pytfa.optim.utils import symbol_sum

from subnetx.io.json import load_json_model, save_json_model

from subnetx.utils.utils import find_blocked_rxns
from subnetx.core.ranking import add_integer_cuts

from addSubnetEcoli import find_target_id
from eval_netws import target_list
from enumerate_subnets import path_load, path_mod, TARGET_RXN_ID
from find_minimal_subnets import product_threshold
from extract_min_pathways import path_save

# User input
new_path_save = '{}/pthws_ids_{}_min.csv'

if __name__ == "__main__":
    
    for target in target_list:
        ### preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        model = load_json_model(path_load.format(target, target))
        data = pd.read_csv(path_save.format(the_path, product_threshold))
        
        pthw_list = []
        for x in data['reactions']:
            rxn_ids = x.replace('[','').replace(']','').replace("'",'').split(', ')
            pthw_list += [rxn_ids]
        
        result = {}
        for ind, rxn_list in enumerate(pthw_list):
            result[ind] = {'reactions': ' // '.join([model.rxn_lexicon[r]['LCSBID'] for r in rxn_list])}
        
        pd.DataFrame.from_dict(result,
                               orient='index').to_csv(
                                   new_path_save.format(the_path,
                                                    product_threshold))
        