#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os.path import join as pjoin
import pandas as pd

import cobra
from copy import deepcopy

from pytfa.io.json import save_json_model, load_json_model

from subnetx.core.chassis import ChassisModel       
from subnetx.io.parser import input_parser_pthw, import_bch_data

from subnetx.core.integration import integrate_pathway

tol = 1e-4 # a tolerance for the product yield

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
solver = GUROBI

target2hub = pd.read_excel('../data/Compound-precursor.xlsx',
                             sheet_name='Sheet1',
                             header=0)

def load_model(host):
    if host == 'ecoli':
        cobra_model = load_json_model('../models/iJO1366.json')
        cobra_model.reactions.EX_glc__D_e.lower_bound = -10
        cobra_model.reactions.BIOMASS_Ec_iJO1366_WT_53p95M.lower_bound = 0
        cobra_model.reactions.EX_co2_e.lower_bound = 0 # blocking Co2 uptake
    elif host == 'yeast':
        cobra_model = load_json_model('../models/Yeast8.json')
        cobra_model.reactions.r_1714.lower_bound = -10
        cobra_model.reactions.r_4041.lower_bound = 0
        cobra_model.reactions.r_1672.lower_bound = 0 # blocking Co2 uptake
    return cobra_model

def load_tmodel(host):
    if host == 'ecoli':
        tmodel = load_json_model('../models/ecoliTFA.json')
        tmodel.reactions.EX_glc__D_e.lower_bound = -10
        tmodel.reactions.BIOMASS_Ec_iJO1366_WT_53p95M.lower_bound = 0
        tmodel.reactions.EX_co2_e.lower_bound = 0 # blocking Co2 uptake
    elif host == 'yeast':
        tmodel = load_json_model('../models/yeastTFA.json')
        tmodel.reactions.r_1714.lower_bound = -10
        tmodel.reactions.r_4041.lower_bound = 0
        tmodel.reactions.r_1672.lower_bound = 0 # blocking Co2 uptake
    tmodel.solver.problem.Params.FeasibilityTol = 1e-9 # it is MILP and solver accuracy should be higher 
    return tmodel

def load_met_lexicon(host):
    if host == 'ecoli':
        lexicon = pd.read_csv('../data/LCSBids_iJO1366.csv', index_col = 0)
    elif host == 'yeast':
        lexicon = pd.read_csv('../data/LCSBids_Yeast8.csv', index_col = 0)
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

# main function:
def model_builder(the_path, target, hosts = ['ecoli','yeast'],
                  branched = False, thermo_eval = True):
    # the_path: ~/Desktop/GITO/patheval/oldForm/results_processed-master/EPFL/Target
    
    # the_path is the directory where the data for this compound exist
    ## this part is only target specific.
    # this function extract data for the linear pathways
    met_list, rxn_list, pthw_list = input_parser_pthw(the_path)
    bchd_pthw_list = import_bch_data(the_path) if branched else None   
    
    for org in hosts:
        sol_list = []
        
        for ind, this_pthw in pthw_list.iterrows():
            ### this part is host specific in addition.
            raw_model = load_model(org)
            met_lexicon = load_met_lexicon(org).to_dict(orient='index')
            # rxn_lexicon = load_rxn_lexicon(org).to_dict(orient='index')
            cobra_model = ChassisModel(raw_model, organism=org, met_lexicon=met_lexicon,
                                    rxn_lexicon={}, inplace=True)

            # If the precursor is not in the model, the supporting pathways should be added
            which_hubs = find_which_hub(target)
            for hub_ind, hub in enumerate(which_hubs):
                hub_path = the_path.replace(target, hub)
                hub_met_list, hub_rxn_list, hub_pthw_list = input_parser_pthw(hub_path)
                hub_bchd_pthw_list = import_bch_data(hub_path) if branched else None
                met_prefix = 'MET_{}_'.format(hub_ind)
                rxn_prefix = 'RXN_{}_'.format(hub_ind)
                integrate_pathway(cobra_model, hub_met_list, hub_rxn_list,hub_pthw_list, hub_bchd_pthw_list,
                            met_id_prefix = met_prefix, rxn_id_prefix = rxn_prefix)
            
            ### this part is pathway specific in addition.
            target_rxn = integrate_pathway(cobra_model, met_list, rxn_list,
                                     this_pthw, bchd_pthw_list)
            cobra_model.objective = target_rxn
            # check if with FBA, the product can be produced
            fba_sol = cobra_model.slim_optimize()
            # file_path = the_path + '/' + cobra_model.name + pthw_id + '.json'
            # save_json_model(cobra_model, file_path)
            if fba_sol < tol or not thermo_eval:
                tfa_sol = 0
            else:
                # load the thermo model instead of constructing it
                raw_tmodel = load_tmodel(org)
                tfa_model = ChassisModel(raw_tmodel, organism=org, met_lexicon=met_lexicon,
                                    rxn_lexicon={}, inplace=True)
                # If the precursor is not in the model, the supporting pathways should be added
                for hub_ind, hub in enumerate(which_hubs):
                    hub_path = the_path.replace(target, hub)
                    hub_met_list, hub_rxn_list, hub_pthw_list = input_parser_pthw(hub_path)
                    met_prefix = 'MET_{}_'.format(hub_ind)
                    rxn_prefix = 'RXN_{}_'.format(hub_ind)
                    integrate_pathway(tfa_model, hub_met_list, hub_rxn_list,hub_pthw_list, hub_bchd_pthw_list,
                            met_id_prefix = met_prefix, rxn_id_prefix = rxn_prefix)
                    
                target_rxn = integrate_pathway(tfa_model, met_list, rxn_list,
                                     this_pthw, bchd_pthw_list)
                tfa_model.prepare()
                tfa_model.convert()
                # file_path = the_path + '/' + tfa_model.name + pthw_id + '.json'
                tfa_model.solver = solver
                tfa_model.objective = target_rxn
                tfa_sol = tfa_model.slim_optimize()
                # save_json_model(tfa_model, file_path)
                
            sol_list += [(this_pthw['P_PR_UID'], fba_sol, tfa_sol)]
            data = pd.DataFrame(sol_list, columns = ['P_PR_UID', 'FBA', 'TFA'])
            data.to_csv('{}/pathway_eval_{}.csv'.format(the_path,org))
