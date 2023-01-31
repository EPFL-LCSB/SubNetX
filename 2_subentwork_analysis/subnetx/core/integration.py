#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 15:26:11 2020

@author: omid
"""
import pandas as pd

from ..io.parser import pthw_parser, merge_pthw_info, find_pthw

RXN_ID_COL = 'R_PR_UID'

def integrate_pathway(model, met_list, rxn_list, pthw, bchd_pthw_list,
                 met_id_prefix = 'metabolite_', rxn_id_prefix = 'reaction_'):
    '''
    This function adds a pathway (linear or branched) or a set of parallel pathways

    Parameters
    ----------
    model : HostModel
        The model for the host.
    met_list : pd.DataFrame
        The information for all the metabolites that should be added.
    rxn_list : pd.DataFrame
        The information for all the reactions that should be added.
    pthw : pd.Series or pd.DataFrame
        The information of the pathway(s).

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
        
    # extract the data from the files in a formatted way
    if isinstance(pthw, pd.Series): # a single (linear) pathway
       rxns, target, boundary_mets, pthw_id, bchd_pthws = pthw_parser(pthw)
    elif isinstance(pthw, pd.DataFrame): # a set of (branched) pathways
       rxns, target, boundary_mets, pthw_id, bchd_pthws = merge_pthw_info(pthw)
       
    # adding supporting pathways to synthesize cosubstrates
    while len(bchd_pthws) > 0: # it is possible that branched pathways have also branched pathways!
       pthw_to_be_added = find_pthw(bchd_pthws, bchd_pthw_list)
       bchd_rxns, _, bchd_boundary_mets, bchd_pthw_id, bchd_pthws = \
           merge_pthw_info(pthw_to_be_added)
       rxns = list(set(rxns + bchd_rxns))
       boundary_mets = list(set(boundary_mets + bchd_boundary_mets))
       pthw_id = list(set(pthw_id + bchd_pthw_id))
    
    
    # adding the reactions and metabolites for the pathway
    model.add_mets_rxns(met_list, rxn_list, rxns, met_id_prefix, rxn_id_prefix) # this also updates lexicons
    
    target_rxn = model.add_target(target) # this also updates the product list
    
    model.add_boundary_rxns(boundary_mets)
    
    return target_rxn


def integrate_network(model, met_list, rxn_list, target,
                 met_id_prefix = 'metabolite_', rxn_id_prefix = 'reaction_'):
    '''
    This function adds a pathway (linear or branched) or a set of parallel pathways

    Parameters
    ----------
    model : HostModel
        The model for the host.
    met_list : pd.DataFrame
        The information for all the metabolites that should be added.
    rxn_list : pd.DataFrame
        The information for all the reactions that should be added.
    target : str
        The (LCSB) ID for the target metabolites.
    

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    rxns = rxn_list[RXN_ID_COL].to_list()
        
    # adding the reactions and metabolites for the pathway
    model.add_mets_rxns(met_list, rxn_list, rxns, met_id_prefix, rxn_id_prefix) # this also updates lexicons
    
    target_rxn = model.add_target(target) # this also updates the product list
        
    return target_rxn