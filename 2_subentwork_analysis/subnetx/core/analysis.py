#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 19:02:30 2021

@author: omid
"""
from warnings import warn


def integrate_strategy(raw_model, integrated_model, rxn_list, target_rxn):
    '''
    A function to integrate specific reactions for a strategy, pathway or subnet.

    Parameters
    ----------
    raw_model : ChassisModel
    integrated_model : ChassisModel with the integrated whole subnetwork
    rxn_list : list of rxn ids to be added to the raw model
    target_rxn : str

    Returns
    -------
    None.

    '''
    if target_rxn not in rxn_list:
        warn('The target reaction was not in the reaction list to be added.')
        rxn_list.append(target_rxn)
    
    for rxn_id in rxn_list:
        rxn_add = integrated_model.reactions.get_by_id(rxn_id)
        raw_model.add_reaction(rxn_add)
        
    
    
    # updating the lexicons
    raw_model.met_lexicon = {met.id:integrated_model.met_lexicon[met.id] for \
                             met in raw_model.metabolites}
    raw_model.rxn_lexicon = {rxn.id:integrated_model.rxn_lexicon[rxn.id] for \
                             rxn in raw_model.reactions if rxn.id != target_rxn} # the demand for the target does not have LCSBID
    raw_model.objective = target_rxn
    
    return