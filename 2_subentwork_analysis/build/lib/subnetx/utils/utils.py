#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 11:39:03 2020

@author: omid
"""

import cobra
import pandas as pd
import numpy as np
import requests
from warnings import warn


kegg_url = 'http://rest.kegg.jp/get/{}/'

def get_molecular_weight(kegg_id):
    '''
    

    Parameters
    ----------
    kegg_id : str

    Returns
    -------
    mol_wieght : float
        the molecular weight in g/mol.

    '''
    response = requests.post(kegg_url.format(kegg_id))
    
    if response.ok:
        eol_ix = response.text.split('\n')
        try:
            mol_wieght_field = [x for x in eol_ix if 'MOL_WEIGHT' in x][0]
            _, mol_wieght = mol_wieght_field.split('  ')
            mol_wieght = float(mol_wieght)
        except IndexError:
            warn('This target compound does not have an assigned kegg ID and \
              molecular weight should be found by another means')
            mol_wieght = 100
        
    else:
        mol_wieght = np.nan
    
    return mol_wieght

def find_blocked_rxns(model, check_rxn_id = ''):
    # a function to find the reactions that are  blocked
    
    blocked_rxns = []
    for rxn in model.reactions:
        if check_rxn_id in rxn.id:
            model.objective = rxn
            model.objective_direction = 'min'
            lb = model.slim_optimize()
            model.objective_direction = 'max'
            ub = model.slim_optimize()
            if ub == 0 and lb==0:
                blocked_rxns.append(rxn)
    
    model.remove_reactions(blocked_rxns)
    return blocked_rxns