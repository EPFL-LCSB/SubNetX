#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 19:06:12 2021

@author: omid
"""
import os
import numpy as np
import pandas as pd

from tqdm import tqdm

from pytfa.optim.variables import ForwardBackwardUseVariable
from pytfa.optim.utils import symbol_sum

from subnetx.io.json import load_json_model, save_json_model

from subnetx.core.ranking import add_integer_cuts

from add_subnet import find_target_id
from eval_netws import target_list
from enumerate_subnets import path_save, path_mod, TARGET_RXN_ID

# User input
product_threshold = 0.5
num_pthw = 10


path_load = path_save
path_enum_pthws = path_mod + '/{}/minimal_subnets_{}_production.csv'

            
def enumerator(index, model, target_rxn):
    
    model.objective_direction = 'min'
    model.slim_optimize()
    
    # Also taking the values of the binary variables
    b_use_vars = model.get_variables_of_type(ForwardBackwardUseVariable)
    
    try:
        ret = {
        # 'index' : index,
        'available_substrate' : - model.reactions.EX_glc__D_e.flux,
        }
        ret['product_flux'] = model.reactions.get_by_id(target_rxn).flux
        
        for b_var in b_use_vars:
            ret[b_var.name] = b_var.variable.primal
                
        # Adding the integer cut constraints to enumerate the other solutions
        add_integer_cuts(model, b_use_vars, index)
        
    except Exception: # runtime error happened -> no solution is found in this time 
        model.solver.configuration.timeout = 1 # change the time to avoid wasting more time
        ret = {
        'available_substrate' : np.nan,
        'product_flux'        : np.nan
        }
        for b_var in b_use_vars:
            ret[b_var.name] = np.nan
    
    
    
    sol = pd.Series(ret)
    
    print('The {}th subnet is found.'.format(index))
    return sol


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
            
            
        target_id = find_target_id(target)
        target_met_id = [k for k,v in model.met_lexicon.items() \
                         if v['LCSBID'] == target_id][0]
        target_rxn = TARGET_RXN_ID + target_met_id
        
        ### Find alternative minimal solutions
        model.objective = target_rxn
        max_prod = model.slim_optimize()
        model.reactions.get_by_id(target_rxn).lower_bound = product_threshold * max_prod
        
        obj_vars = [var.variable for var in \
                    model.get_variables_of_type(ForwardBackwardUseVariable)]    
        model.objective = symbol_sum(obj_vars)
        model.objective_direction = 'min'
        
        model.solver.configuration.tolerances.feasibility = 1e-9
        model.solver.configuration.timeout = 4*3600
        
        # Tightening the bounds to increase the solver efficiency
        for rxn in model.reactions:
            if rxn.lower_bound == -1000:
                rxn.lower_bound = -50
            if rxn.upper_bound == 1000:
                rxn.upper_bound = 50
        
            
        pthw_index = pd.Series(np.arange(0, num_pthw, 1))
        data = pthw_index.apply(enumerator , args=[model, target_rxn])
        data.to_csv(path_enum_pthws.format(target, product_threshold))
    