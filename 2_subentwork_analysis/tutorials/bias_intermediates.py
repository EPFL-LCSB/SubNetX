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

from subnetx.utils.utils import find_blocked_rxns
from subnetx.core.ranking import add_integer_cuts, bias_passing_intermediate

from addSubnetEcoli import find_target_id
from eval_netws import target_list
from enumerate_subnets import path_save, path_mod, TARGET_RXN_ID

# User input
product_threshold = 1
num_pthw = 10


path_load = path_save
path_enum_pthws = path_mod + '/{}/biased_subnets_{}_production.csv'

            
def enumerator(index, model, target_rxn):
    
    model.objective_direction = 'min'
    model.slim_optimize()
    
    ret = {
        # 'index' : index,
        'available_substrate' : - model.reactions.EX_glc__D_e.flux,
        }
    
    ret['product_flux'] = model.reactions.get_by_id(target_rxn).flux
    
    # Also taking the values of the binary variables
    b_use_vars = model.get_variables_of_type(ForwardBackwardUseVariable)
    for b_var in b_use_vars:
        ret[b_var.name] = b_var.variable.primal
            
    # Adding the integer cut constraints to enumerate the other solutions
    add_integer_cuts(model, b_use_vars, index)
    
    sol = pd.Series(ret)
    
    print('The {}th subnet is found.'.format(index))
    return sol

def find_targets_in_subnet(model, id_list, metabolites=[]):
    # to find all shiki targets in the subnetwork
    intermediates = []
    if len(metabolites) == 0:
        metabolites = model.metabolites
        
    for met in metabolites:
        if model.met_lexicon[met.id]['LCSBID'] in id_list: # it is a shiki target
            intermediates += [met]
    return intermediates


if __name__ == "__main__":
    
    shiki_100_compounds = pd.read_excel('../data/Compound-precursor.xlsx',
                              sheet_name='Sheet1',
                              header=0)['Compound'].tolist()
    target_lcsb_ids = [find_target_id(x.replace(' ', '_')) for x in shiki_100_compounds]
    
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
        
        # Setting an objective, though it will be modified later
        obj_vars = [var.variable for var in \
                    model.get_variables_of_type(ForwardBackwardUseVariable)]    
        model.objective = symbol_sum(obj_vars)
        model.objective_direction = 'min'
        
        intermediates = find_targets_in_subnet(model, target_lcsb_ids)
        bias_passing_intermediate(model, intermediates, weight = -0.1) # The weight should be negative to reward passing through the intermediate
        
        model.solver.configuration.tolerances.feasibility = 1e-9
        model.solver.configuration.timeout = 2*3600
        pthw_index = pd.Series(np.arange(0, num_pthw, 1))
        data = pthw_index.apply(enumerator , args=[model, target_rxn])
        data.to_csv(path_enum_pthws.format(target, product_threshold))
    