#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 19:02:30 2021

@author: omid
"""
import pickle
import zlib

from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
from pytfa import ThermoModel

from warnings import warn
import os
from ..thermo.thermodb import append_thermodbs

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

def thermo_evaluation(model, thermo_data=None, thermo_data_path=None,
                      new_metabolites=False, new_met_list=[]):
    
    # Preparing the thermo-database
    if thermo_data is None:
        if thermo_data_path is None or \
        not os.path.exists(thermo_data_path): # online thermo database
            warn('The thermo database is not provided so eQuilibrator is runnning.')
            # do this only for non-new metabolites as it saves for later use
            met_list = [met for met in model.metabolites if met not in new_met_list]
            thermo_data = build_thermo_from_equilibrator(met_list)
            # in case the path is provided, try to save the thermo database for future use
            if thermo_data_path is not None:
                warn('The thermo database is being saved in the reference provided.')
                with open(thermo_data_path, 'wb') as file:
                    # Convert it to string with Pickle
                    data = pickle.dumps(thermo_data, pickle.HIGHEST_PROTOCOL)
                    # Compress it with zlib
                    data = zlib.compress(data)
                    # And write it to disk
                    file.write(data)
                
        else: # the reference to load a local thermodatabase
            from pytfa.io import load_thermoDB
            thermo_data = load_thermoDB(thermo_data_path)
            
            
    if new_metabolites and new_met_list: # even in the case of local database we may still have new metabolites
        thermo_data_subset = build_thermo_from_equilibrator(new_met_list)
        thermo_data = append_thermodbs(thermo_data, thermo_data_subset)
    
    # Creating the thermo model
    tmodel = ThermoModel(thermo_data, model)
    # since the next functions try to match seed_id to the thermo database
    for met in tmodel.metabolites:
        met.annotation.update({'seed_id' : met.id})
        
    tmodel.prepare()
    tmodel.convert()           
    
    return thermo_data, tmodel