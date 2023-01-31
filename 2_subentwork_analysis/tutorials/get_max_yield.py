#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:31:17 2020

@author: omid
"""
from os.path import join as pjoin
import pandas as pd
import os.path
import numpy as np
from warnings import warn

from thermoModelbuilder import load_model
from subnetx.utils.utils import get_molecular_weight


path_mod = '../work/results_processed-master/EPFL'
# Compound data
target_list = pd.read_excel('../data/Compound-precursor.xlsx',
                             sheet_name='Sheet1',
                             header=0)['Compound'].tolist()
# Host data
host_list = ['ecoli','yeast']
host_glc_id = {'ecoli' : 'EX_glc__D_e',
               'yeast' : 'r_1714'}
host_model_dict = {org:load_model(org) for org in host_list}

host_glc_uptake = {org:model.reactions.get_by_id(host_glc_id[org]).lower_bound \
                   for org, model in host_model_dict.items()}

GLC_MW = 180.156 # glucose molecular weight g/mol
GLC_FLUX = 10 # mmol/(gDW.h)

def get_target_mol_weight(the_path):
    met_list = pd.read_csv(pjoin(the_path,'metabolites.tsv'), sep='\t')
    pthw_list = pd.read_csv(pjoin(the_path,'pathways.tsv'), sep='\t')
    kegg_id = \
        met_list[met_list['M_PR_UID']==pthw_list['P_PR_TARGET'][0]]['M_XR_KEGG'].values[0]
    if pd.isna(kegg_id):
        warn('This target compound does not have an assigned kegg ID and \
              molecular weight should be found by another means')
        return 100
    mol_weight = get_molecular_weight(kegg_id)
    return mol_weight
    
if __name__ == "__main__":
    modifiers = ['max FBA yield {}', 'max TFA yield {}', 
                 'pct FBA feas. pathways {}', 'pct TFA feas. pathways {}']
    columns = [x.format(y)for y in host_list for x in modifiers]
    result = pd.DataFrame([],columns = columns, index = [t.replace(' ', '_') \
                                                         for t in target_list])
    for target in target_list:
        # preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            continue
        # load flux data
        result_dict = {host:pd.read_csv('{}/pathway_eval_{}.csv'.format(the_path,host)
                                        ,header=0,index_col=0) \
                       for host in host_list}
        for org, outflux in result_dict.items():
            glc_flux = GLC_FLUX # abs(host_glc_uptake[org]) # positive flux
            target_mol_weight = get_target_mol_weight(the_path)
            # the yield is: (target_flux*target_weight)/(glc_flux*glc_weight)
            FBA_yield = outflux['FBA'] * target_mol_weight/\
                (glc_flux * GLC_MW)
            TFA_yield = outflux['TFA'] * target_mol_weight/\
                (glc_flux * GLC_MW)
            yields = pd.concat([outflux['P_PR_UID'],FBA_yield,TFA_yield],axis=1)
            yields.to_csv('{}/pathway_yield_{}.csv'.format(the_path,org))
            # for a general overview, calculate max yield and 
            # the percent of feasible pathways for each compound 
            fba_list = FBA_yield.tolist()
            tfa_list = TFA_yield.tolist()
            
            result[modifiers[0].format(org)][target] = max(fba_list)
            result[modifiers[1].format(org)][target] = max(tfa_list)
            result[modifiers[2].format(org)][target] = np.count_nonzero(fba_list)/len(fba_list)
            result[modifiers[3].format(org)][target] = np.count_nonzero(tfa_list)/len(tfa_list)
    
    result.to_csv('Overview.csv')
            