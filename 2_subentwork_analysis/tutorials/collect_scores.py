#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 17:16:44 2021

@author: omid
"""
import os
import pandas as pd
from tqdm import tqdm
import pytfa
from pytfa.io import load_thermoDB
from subnetx.core.chassis import ChassisModel, LCSBID
from subnetx.core.ranking import assign_bin_var2rxn
from subnetx.core.analysis import integrate_strategy
from subnetx.io.json import load_json_model
from addSubnetEcoli import load_model, load_cobra_model, find_target_id
from eval_netws import target_list, path_mod
from enumerate_subnets import TARGET_RXN_ID, HETERO_RXN_ID, PRE_RXN_ID
from enumerate_subnets import path_save as integ_model_path
from extract_min_pathways import path_save as path_minimal
# from bias_intermediates import find_targets_in_subnet
from find_high_score import bridgit_score
from extract_biased_pathways import path_save as path_biased

minimal_prod_thresh = [0.25,0.5,0.75,1]
hscore_prod_thresh = [0.25,0.5,0.75,1] 
thermodb_path = 'results_processed-master/EPFL/thermo_all/thermodata_ecoli.thermodb'


def prep_score_table(target, rxn_set, pathway_list):
    '''
    

    Parameters
    ----------
    target : str
    rxn_set : list
        the cumulated list of all reactions to be added without redundnacy
    pathway_list : list
        the list of found pathways (each element is a list of reactions) without redundnacy

    Returns
    -------
    ret : pd.DataFrame
        scoring table.

    '''
    # loading data
    
    # loading the models
    raw_model = load_model()
    model = ChassisModel(raw_model, organism='ecoli', met_lexicon={},
                        rxn_lexicon={}, inplace=True)
    
    integ_model = load_json_model(integ_model_path.format(target,target))
    
    # find the production of the product
    target_id = find_target_id(target)
    target_met_id = [k for k,v in integ_model.met_lexicon.items() \
                     if v[LCSBID] == target_id][0]
    target_rxn = TARGET_RXN_ID + target_met_id
    # target reactions is not already in the lexicon
    integ_model.rxn_lexicon[target_rxn] = 'nan'
    
    # Add the reactions altogether        
    integrate_strategy(model, integ_model, rxn_set, target_rxn)
    
    # Then we block all the added reactions
    for rxn in rxn_set:
        model.reactions.get_by_id(rxn).bounds =(0,0)
        
    # open the target reaction so if it is not in the list of reactions of the pathway, we have it
    model.reactions.get_by_id(target_rxn).bounds = (0, 1000)
    
    rxn_dict = {r:integ_model.rxn_lexicon[r] for r in rxn_set}
    
    ret = {}
    # calculate the scores
    for pathway in tqdm(pathway_list, desc='checking pathways'):
        key = [rxn_dict[x][LCSBID] for x in pathway] # a list of rxn LCSB IDs for this pathway
        
        this_dict = {'Size' : len(pathway)} # initiate a dict to store scores for this pathway
        
        pathway_mets = [] # this will be used later for counting the intermediates
        # Unblocking the rxns of this pathway
        for rxn_id in pathway:
            rxn = model.reactions.get_by_id(rxn_id)
            pathway_mets += [met for met in  rxn.metabolites if met not in pathway_mets]
            rxn.bounds = (-1000, 1000)
        
        # finding the yield 
        model.slim_optimize()
        molar_yield = abs(model.objective.value/model.reactions.EX_glc__D_e.flux)
        this_dict['Molar Yield'] = molar_yield 
        this_dict['Gram Yield'] = molar_yield * \
            (model.metabolites.get_by_id(target_met_id).formula_weight / model.metabolites.glc__D_e.formula_weight)
        
        # blocking the rxns of this pathway again
        for rxn_id in pathway:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.bounds = (0, 0)
        
        # counting shiki intermediates this pathway passes through
        total_score = 0
        for rxn_id in pathway:
            lcsb_id  = integ_model.rxn_lexicon[rxn_id][LCSBID] # first find the LSCB ID
            this_score = bridgit_score[lcsb_id] # then bridgit score
            total_score += this_score 
        this_dict['BridgIT score'] = total_score
        
        ret[str(key)] = this_dict
    
    
    # thermo evaluation
    raw_model = load_cobra_model()
    thermodb = load_thermoDB(thermodb_path)
    integrate_strategy(raw_model, integ_model, rxn_set, target_rxn)
    tmodel = pytfa.ThermoModel(thermodb, raw_model)
    tmodel.prepare()
    tmodel.convert()
    
    tmodel = ChassisModel(tmodel, organism='ecoli', met_lexicon=integ_model.met_lexicon,
                        rxn_lexicon=integ_model.rxn_lexicon, inplace=True)
    # Then we block all the added reactions
    for rxn in rxn_set:
        tmodel.reactions.get_by_id(rxn).bounds =(0,0)
        
    # open the target reaction so if it is not in the list of reactions of the pathway, we have it
    tmodel.reactions.get_by_id(target_rxn).bounds = (0, 1000)
    
    for pathway in tqdm(pathway_list, desc='thermodynamic evaluation'):
        key = [rxn_dict[x][LCSBID] for x in pathway] # a list of rxn LCSB IDs for this pathway
        
        # Unblocking the rxns of this pathway
        for rxn_id in pathway:
            rxn = tmodel.reactions.get_by_id(rxn_id)
            rxn.bounds = (-1000, 1000)
            
        thermo_yield = tmodel.slim_optimize()
        ret[str(key)]['Thermodynamic Eval.'] = 'Feasible' if thermo_yield > 0 \
            else 'InFeasible'
        
        # blocking the rxns of this pathway again
        for rxn_id in pathway:
            rxn = tmodel.reactions.get_by_id(rxn_id)
            rxn.bounds = (0, 0)
    
    return ret


if __name__ == "__main__":
    
    for target in target_list:
        target = target.replace(' ', '_')
        the_path = '{}/{}/results'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        path_load = [path_minimal.format(the_path,x) for x in minimal_prod_thresh] \
                    + [path_biased.format(the_path,x) for x in hscore_prod_thresh]
            
        
        # combining the data
        all_data = None 
        for path in path_load:
            data = pd.read_csv(path)
            if all_data is None:
                all_data = pd.DataFrame(columns= data.columns) # initiate it
            all_data = pd.concat([all_data, data])
        
        # constructing a set of reactions to be added from all subntes
        rxn_set = []
        pathway_list = []
        for rxns in all_data['reactions']:
            rxn_list = \
                rxns.replace("'","").replace('[','').replace(']','').split(', ')
            if rxn_list not in pathway_list:
                pathway_list += [rxn_list]
            rxn_set += [x for x in rxn_list if x not in rxn_set]
            
        score_table = prep_score_table(target, rxn_set, pathway_list)    
        table = pd.DataFrame.from_dict(score_table, orient='index')
        table.to_csv('{}/{}/final_scores.csv'.format(path_mod,target))
