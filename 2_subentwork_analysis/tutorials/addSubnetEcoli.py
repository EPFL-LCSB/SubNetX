#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os.path import join as pjoin
import pandas as pd
import cobra

from pytfa.io.json import load_json_model

from subnetx.core.chassis import ChassisModel       
from subnetx.io.parser import input_parser_pthw, input_parser_netw

from subnetx.core.integration import integrate_pathway, integrate_network

from subnetx.io.json import save_json_model

tol = 1e-4 # a tolerance for the product yield

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
solver = GUROBI

ORGANISM = 'ecoli'

target2hub = pd.read_excel('../data/Compound-precursor.xlsx',
                             sheet_name='Sheet1',
                             header=0)
target2lcsb = pd.read_csv('../data/target_LCSBID.csv',
                             header=0)

def load_model():
    cobra_model = load_json_model('../models/iJO1366.json')
    cobra_model.reactions.EX_glc__D_e.lower_bound = -10
    cobra_model.reactions.BIOMASS_Ec_iJO1366_WT_53p95M.lower_bound = 0
    cobra_model.reactions.EX_co2_e.lower_bound = 0 # blocking Co2 uptake
    return cobra_model

def load_tmodel():
    tmodel = load_json_model('../models/ecoliTFA.json')
    tmodel.reactions.EX_glc__D_e.lower_bound = -10
    tmodel.reactions.BIOMASS_Ec_iJO1366_WT_53p95M.lower_bound = 0
    tmodel.reactions.EX_co2_e.lower_bound = 0 # blocking Co2 uptake
    
    tmodel.solver.configuration.tolerances.feasibility = 1e-9 # it is MILP and solver accuracy should be higher 
    return tmodel

def load_cobra_model():
    cobra_model = cobra.io.load_json_model('../models/iJO1366.json')
    cobra_model.reactions.EX_glc__D_e.lower_bound = -10
    cobra_model.reactions.BIOMASS_Ec_iJO1366_WT_53p95M.lower_bound = 0
    cobra_model.reactions.EX_co2_e.lower_bound = 0 # blocking Co2 uptake
    cobra_model.solver.configuration.tolerances.feasibility = 1e-9 # it is MILP and solver accuracy should be higher 
    return cobra_model
    

def preprocess_model(model):
    '''
    A function to correct some mass balances that are found to be problematic.

    '''
    # GLCS1 is a reaction that converts glucose to glycogene: adpglc_c --> adp_c + glycogen_c + h_c
    # It should be corrected as follows: 5.0 adpglc_c --> 5.0 adp_c + glycogen_c + 5.0 h_c
    adpglc_c = model.metabolites.adpglc_c
    adp_c = model.metabolites.adp_c
    h_c = model.metabolites.h_c
    model.reactions.GLCS1.add_metabolites({adpglc_c:-4, adp_c:4, h_c:4})

def load_met_lexicon():
    lexicon = pd.read_csv('../data/LCSBids_iJO1366.csv', index_col = 0)
    return lexicon.astype('Int64').astype(str)

def load_rxn_lexicon():
    lexicon = pd.read_csv('../data/rxn_LCSBID_iJO1366.csv', index_col = 0)
    return lexicon.astype('Int64').astype(str)

def find_which_hub(target):
    which_hubs = list()
    
    this_target2hub = target2hub[target2hub['Compound']==target.replace('_',' ')] 
    while not this_target2hub['Main hub'].hasnans: # precursor is not in the model
        this_hub = this_target2hub['Main hub'].any()
        which_hubs.append(this_hub.replace(' ', '_'))
        this_target2hub = \
            target2hub[target2hub['Compound']==this_hub]
    return which_hubs

def find_target_id(target):
    try:
        tar = target2lcsb[target2lcsb['Compound']==target]['LCSBID'].tolist()[0]
        if isinstance(tar, int):
            tar = str(tar)
        return tar
    except IndexError:
        raise ValueError('The LCSB ID for this coumpound is not automatically found.'
                         'Please add it manually.')

# main function:
def model_builder(the_path, target, thermo_eval = False):
    
    ### Integration
    # the_path: ~/Desktop/GITO/patheval/oldForm/results_processed-master/EPFL/Target
    # the_path is the directory where the data for this compound exist
    
    # this function extract data for the subnetwork
    met_list, rxn_list = input_parser_netw(the_path)
    target_id = find_target_id(target)
    
    sol_list = []
    met_lexicon = load_met_lexicon().to_dict(orient='index')
    rxn_lexicon = load_rxn_lexicon().to_dict(orient='index')
    
    raw_model = load_model()
    cobra_model = ChassisModel(raw_model, organism=ORGANISM, met_lexicon=met_lexicon,
                            rxn_lexicon=rxn_lexicon, inplace=True)

    # If the precursor is not in the model, the supporting pathways should be added
    # which_hubs = find_which_hub(target)
    # for hub_ind, hub in enumerate(which_hubs):
    #     hub_path = the_path.replace(target, hub)
    #     hub_met_list, hub_rxn_list = input_parser_netw(hub_path)
    #     met_prefix = 'MET_{}_'.format(hub_ind)
    #     rxn_prefix = 'RXN_{}_'.format(hub_ind)
    #     integrate_network(cobra_model, hub_met_list, hub_rxn_list, find_target_id(hub),
    #                 met_id_prefix = met_prefix, rxn_id_prefix = rxn_prefix)
    
    ## integrating the network
    target_rxn = integrate_network(cobra_model, met_list, rxn_list, target_id)
    cobra_model.objective = target_rxn
    
    ### Stoichiometric evaluation
    # solve the issue with the unbalanced reactions
    preprocess_model(cobra_model)
    # check if with FBA, the product can be produced
    fba_sol = cobra_model.slim_optimize()
    
    
    ### Saving the model
    file_path = the_path + '/' + cobra_model.name + '_' + target + '.json'
    save_json_model(cobra_model, file_path)
    
    ### Thermo evaluation
    if fba_sol < tol or not thermo_eval:
        tfa_sol = 0
    else:
        # load the thermo model instead of constructing it
        raw_tmodel = load_tmodel()
        tfa_model = ChassisModel(raw_tmodel, organism=ORGANISM, met_lexicon=met_lexicon,
                            rxn_lexicon=rxn_lexicon, inplace=True)
        # If the precursor is not in the model, the supporting pathways should be added
        # for hub_ind, hub in enumerate(which_hubs):
        #     hub_path = the_path.replace(target, hub)
        #     hub_met_list, hub_rxn_list = input_parser_netw(hub_path)
        #     met_prefix = 'MET_{}_'.format(hub_ind)
        #     rxn_prefix = 'RXN_{}_'.format(hub_ind)
        #     integrate_network(tfa_model, hub_met_list, hub_rxn_list, find_target_id(hub),
        #             met_id_prefix = met_prefix, rxn_id_prefix = rxn_prefix)
        
        ## integrating the network
        target_rxn = integrate_network(tfa_model, met_list, rxn_list, target_id)
        
        tfa_model.prepare()
        tfa_model.convert()
        tfa_model.solver = solver
        tfa_model.objective = target_rxn
        tfa_sol = tfa_model.slim_optimize()
        # file_path = the_path + '/' + tfa_model.name + pthw_id + '.json'
        # save_json_model(tfa_model, file_path)
    
    ### Saving the solution    
    sol_list += [(fba_sol, tfa_sol)]
    data = pd.DataFrame(sol_list, columns = [ 'FBA', 'TFA'])
    data.to_csv('{}/network_eval_{}.csv'.format(the_path,ORGANISM))
    return cobra_model
