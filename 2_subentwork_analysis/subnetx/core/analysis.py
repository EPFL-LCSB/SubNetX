#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 19:02:30 2021

@author: omid
"""
import pickle
import zlib
import pandas as pd
import ast
from openbabel import openbabel

from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
from pytfa import ThermoModel
from pytfa.io import load_thermoDB
from pytfa.utils.numerics import BIGM_DG

from warnings import warn
import os
from ..thermo.thermodb import append_thermodbs, build_thermo_novel

import secrets
import string

import sys
import warnings
from contextlib import contextmanager

# Get the directory where this Python file is located
base_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the full path to the data file
file_path = os.path.join(base_dir, '..', '..', 'data', 'thermo_data.thermodb')

@contextmanager
def suppress_output():
    # Redirect stdout, stderr and suppress warnings
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                yield
            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr


def rand_id(length=6):
    # generate randon hashes
    alphabet = string.ascii_letters + string.digits  # a-zA-Z0-9
    return ''.join(secrets.choice(alphabet) for _ in range(length))

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
    if isinstance(thermo_data, dict): # directly provided
        pass
    
    elif thermo_data is None:
        if thermo_data_path is None or \
        not os.path.exists(thermo_data_path): # online thermo database
            warn('The thermo database is not provided so eQuilibrator is runnning.')
            # do this only for non-new metabolites as it saves for later use
            met_list = [met for met in model.metabolites if met not in new_met_list]
            with suppress_output():
                thermo_data = build_thermo_from_equilibrator(met_list)
            # for these compounds we may also use pytfa thermodb
            # for cases with too large error (usually inorganics), use the other database to replace
            pytfa_thermodb = load_thermoDB(file_path)
            removable_ids = []
            calorie2joule = 4.184 
            metData = thermo_data['metabolites']
            for id_, data in metData.items():
                if data['deltaGf_err'] > BIGM_DG: # too large errors
                    # the other database is keyed by SEED IDs
                    try:
                        seed_id = model.metabolites.get_by_id(id_).annotation['seed.compound']
                        # the SEED IDs in annotaions might be a string, a list or, like "['cpd00009', 'cpd27787']"
                        if isinstance(seed_id, str) and '[' in seed_id:
                            seed_id = ast.literal_eval(seed_id)
                        if isinstance(seed_id, list):
                            seed_id = next(iter(seed_id))
                        new_data = pytfa_thermodb['metabolites'][seed_id]
                        # rare cases the pytfa database includes junk data
                        if (new_data['deltaGf_std'] == 0 and \
                            new_data['deltaGf_err'] == BIGM_DG):
                            removable_ids += [id_]
                        metData[id_]['deltaGf_err'] = \
                            new_data['deltaGf_err'] * calorie2joule
                        metData[id_]['deltaGf_std'] = \
                            new_data['deltaGf_std'] * calorie2joule
                        metData[id_]['pKa'] = new_data['pKa']
                        metData[id_]['formula'] = new_data['formula']    
                    except KeyError:
                        pass
                # correct the charges if needed
                met_charge = model.metabolites.get_by_id(id_).charge
                if data['charge_std'] != met_charge:
                    if data['charge_std'] == 0:
                        metData[id_]['charge_std'] = met_charge
                    else:
                        warn(f"The charge for {id_} may not be properly calculated.")
            thermo_data['metabolites'] = {k:v for k,v in thermo_data['metabolites'].items() \
                                          if k not in removable_ids}
            
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
            thermo_data = load_thermoDB(thermo_data_path)
            
            
    if new_metabolites and new_met_list: # even in the case of local database we may still have new metabolites
        # new metabolites are those that are probably not even cached database
        # we need to create a DataFrame like the following:
        # compound_df = pd.DataFrame(
        #     data=[
        #         ["OC(=O)C1=CC(NC(=O)C2=CC=CC=C2)=C(O)C=C1", "3B4HA", "3-Benzamido-4-hydroxybenzoic acid"],
        #         ["NC1=C(O)C=CC(=C1)C(O)=O", "3A4HA", "3-Amino-4-hydroxybenzoic acid"]
        #     ],
        #     columns=["struct","coco_id", "name"]
        # )
        novel_met_list = [met for met in new_met_list \
                          if met.id not in thermo_data['metabolites']] # only those without thermo calculated
        compound_df = pd.DataFrame(data=[
            [met.annotation['smiles'], rand_id(), met.name] for met in novel_met_list
            ], columns=["struct","coco_id", "name"])
        
        met_smiles_dict = {met.id:met.annotation['smiles'] for met in novel_met_list}
        with suppress_output():
            openbabel.obErrorLog.StopLogging() # Disable the warnings
            thermo_data_subset = build_thermo_novel(compound_df, met_smiles_dict)
        for id_, data in thermo_data_subset['metabolites'].items():
            if data['deltaGf_err'] > BIGM_DG: # too large errors
                data['deltaGf_err'] = 0.2 * data['deltaGf_std'] # if it's too large set the 20%
        thermo_data = append_thermodbs(thermo_data, thermo_data_subset)

    # Creating the thermo model
    if not isinstance(model, ThermoModel):
        tmodel = ThermoModel(thermo_data, model)
    else:
        model.thermo_data = thermo_data # just update the thermo data
        
    # Because the next functions try to match seed_id to the thermo database
    for met in tmodel.metabolites:
        met.annotation.update({'seed_id' : met.id})   
    tmodel.prepare()
    tmodel.convert()           
    
    return thermo_data, tmodel