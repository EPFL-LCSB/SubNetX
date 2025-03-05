#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 19:06:12 2021

@author: omid
"""
import os
import numpy as np
import pandas as pd

from pytfa.optim.variables import ForwardBackwardUseVariable

from subnetx.io.json import load_json_model, save_json_model

from subnetx.utils.utils import find_blocked_rxns
from subnetx.core.ranking import add_integer_cuts, assign_bin_var2rxn

from add_subnet import find_target_id
from eval_netws import target_list, path_mod

# User input
num_pthw = 1
prepared_model = False # if the pruned model with binary variables is saved


path_load = path_mod + '/{}/Model_{}.json'
path_save = path_mod + '/{}/Model_post-analysis_{}.json'
path_enum_pthws = path_mod + '/{}/enumerated_subnets.csv'

HETERO_RXN_ID = 'reaction_'
PRE_RXN_ID = 'RXN_'
TARGET_RXN_ID = 'DM_'

            
def enumerator(index, model, target_rxn):
    
    ret = {
        # 'index' : index,
        'available_substrate' : - model.reactions.r_1714.lower_bound,
        }
    
    model.objective = target_rxn
    model.objective_direction = 'max'
    model.slim_optimize()
    
    ret['product_flux'] = model.objective.value
    
    # Also taking the values of the binary variables
    b_use_vars = model.get_variables_of_type(ForwardBackwardUseVariable)
    for b_var in b_use_vars:
        ret[b_var.name] = b_var.variable.primal
            
    # Adding the integer cut constraints to enumerate the other solutions
    add_integer_cuts(model, b_use_vars, index)
    
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
        if prepared_model:
            model = load_json_model(path_save.format(target, target))
        else:
            model = load_json_model(path_load.format(target, target))
            ### First round of pruning: to remove the reactions that cannot carry flux
            blocked_rxns = find_blocked_rxns(model, check_rxn_id= HETERO_RXN_ID) # this also removes the blocked
            assign_bin_var2rxn(model, HETERO_RXN_ID)
            assign_bin_var2rxn(model, PRE_RXN_ID)
            save_json_model(model, path_save.format(target, target))
            
        target_id = find_target_id(target)
        target_met_id = [k for k,v in model.met_lexicon.items() \
                         if v['LCSBID'] == target_id][0]
        target_rxn = TARGET_RXN_ID + target_met_id
        
        ### Find alternative solution
        model.solver.configuration.tolerances.feasibility = 1e-9
        pthw_index = pd.Series(np.arange(0, num_pthw, 1))
        data = pthw_index.apply(enumerator , args=[model, target_rxn])
        data.to_csv(path_enum_pthws.format(target))
    