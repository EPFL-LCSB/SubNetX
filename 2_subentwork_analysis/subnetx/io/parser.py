#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 20:25:53 2021

@author: omid
"""
import pandas as pd
from os.path import join as pjoin
import os
import re

RXN_ID_COL   = 'R_PR_UID'
RXN_EQ_COL   = 'R_PR_STOICH'
MET_ID_COL   = 'M_PR_UID'
MET_FOR_COL  = 'M_PR_FORMULA'
MET_NAM_COL  = 'M_PR_NAME'
MET_BAR_COL  = 'M_PR_CHARGE'
MET_KEGG_COL = 'M_XR_KEGG'
MET_INCH_COL  = 'M_XR_INCHIKEY' 
MET_SMI_COL  = 'M_XR_SMILES' 

# Get the directory where this Python file is located
base_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the full path to the data file
filename = os.path.join(base_dir, '..', '..', 'data', 'KEGG2SEED_update.xlsx')

def input_parser_netw(the_path):
    # a function to parse the pathways, metabolites and reactions input files
    # target: str, the name of taget compound
    
    met_list = pd.read_csv(pjoin(the_path,'output_optimization_input/compounds.tsv'), sep='\t')
    rxn_list = pd.read_csv(pjoin(the_path,'output_optimization_input/reactions.tsv'), sep='\t')
    # the metabolite IDs are numbers, we should convert them to strings
    met_list['M_PR_UID'] = met_list['M_PR_UID'].astype(str)
    rxn_list['R_PR_UID'] = rxn_list['R_PR_UID'].astype(str)
    
    return met_list, rxn_list

def input_parser_pthw(the_path):
    # a function to parse the pathways, metabolites and reactions input files
    # target: str, the name of taget compound
    
    met_list = pd.read_csv(pjoin(the_path,'output_optimization_input/compounds.tsv'), sep='\t')
    rxn_list = pd.read_csv(pjoin(the_path,'output_optimization_input/reactions.tsv'), sep='\t')
    pthw_list = pd.read_csv(pjoin(the_path,'output_optimization_input/pathways.tsv'), sep='\t')
    # the metabolite IDs are numbers, we should convert them to strings
    met_list['M_PR_UID'] = met_list['M_PR_UID'].astype(str)
    rxn_list['R_PR_UID'] = rxn_list['R_PR_UID'].astype(str)
    
    return met_list, rxn_list, pthw_list

def pthw_parser(pthw):
    r_1 = ' | '
    r_2 = '|'
    rxns = safe_split(pthw['P_PR_REACTIONS'], r_1, r_2) # list
    target = str(pthw['P_PR_TARGET']) # hash
    pthw_id = [pthw['P_PR_UID']] # list
    boundary_mets = [str(x) for x in safe_split(pthw['P_PR_BOUNDARY'], r_1, r_2)] # list
    try:
        bchd_pthws = [x for x in safe_split(pthw['P_PR_BOUNDARY_ORIGIN'], r_1, r_2) \
                       if x != 'Model']
    except KeyError: # The pathway is linear
        bchd_pthws = []
    return rxns, target, boundary_mets, pthw_id, bchd_pthws

# Helper function to parse one side
def parse_side(side, sense):
    # if sense 1 produced if -1 consumed in the reaction
    compounds = side.split('+')
    parsed = {}
    for compound in compounds:
        compound = compound.strip()
        # Match formats: (2)A, 2 A, 1B, A
        match = re.match(r'(?:\((\d+)\)\s*|(\d+)\s+)?(.+)', compound)
        if match:
            coef1, coef2, name = match.groups()
            coef = int(coef1 or coef2 or 1)
            name = name.strip()
            parsed[name] = sense * coef
    return parsed


def rxn_parser(rxns, rxn_list):
    # it must contain the IDs and formulae of reactions
    mets = [] # list
    rxns_stoich = dict() # dict of dict
    remove_candid = []
    for rxn in rxns:
        rxn_dict = dict()
        rxn_data = rxn_list[rxn_list[RXN_ID_COL]==rxn]
        rxn_print = rxn_data[RXN_EQ_COL].values[0]
        participants = safe_split(rxn_print, ' ==> ', ' <=> ', ' <==>') # splitting reactants and products
        reactants = participants[0]
        products = participants[1]
        reactant_dict = parse_side(reactants, -1)
        product_dict = parse_side(products, 1)
        # check if at least a compound is both among products and reactants, not to add the rxn
        if len(set(reactant_dict.keys()) - set(product_dict.keys())) != \
            len(set(reactant_dict.keys())):
                remove_candid.append(rxn)
                continue
            
        for the_met, the_coef in reactant_dict.items():
            if the_met not in mets:
                mets.append(the_met)
            rxn_dict[the_met] = the_coef
        for the_met, the_coef in product_dict.items():
            if the_met not in mets:
                mets.append(the_met)
            rxn_dict[the_met] = the_coef
        rxns_stoich[rxn] = rxn_dict
        
    # removing the spoiled reactions
    for rxn in remove_candid:
        rxns.remove(rxn)
    return mets, rxns_stoich

def met_parser(mets, met_list, host):
    # at least it must contain the IDs and names of metabolites
    mets_formula = dict() # dict
    mets_name = dict() # dict
    mets_charge = dict() # dict
    mets_annotation = dict() # dict
    numerator = 0 # to enumerate metabolites without annotations
    for met in mets:
        met_data = met_list[met_list[MET_ID_COL]==met]
        mets_name[met] = met_data[MET_NAM_COL].values[0] if not pd.isna(met_data[MET_NAM_COL].values[0]) \
            else 'generic_metabolite'
        
        if MET_FOR_COL in met_data.columns:
            mets_formula[met] = met_data[MET_FOR_COL].values[0] \
                if not pd.isna(met_data[MET_FOR_COL].values[0]) else ''
        else:
            mets_formula[met] = ''
            
        if MET_BAR_COL in met_data.columns:
            mets_charge[met] = met_data[MET_BAR_COL].values[0] \
                if not pd.isna(met_data[MET_BAR_COL].values[0]) else 0
        else:
            mets_charge[met] = 0
                
        # Use annotations if exist to improve integration
        if MET_KEGG_COL in met_data.columns: # Specific for KEGG annotation
            with open(filename, 'rb') as kegg2seed:
                kegg_id = met_data[MET_KEGG_COL].values[0]
                if not pd.isna(kegg_id):
                    try:
                        mets_annotation[met] = { 'kegg.compound' : kegg_id, 
                                                'seed.compound' : \
                            kegg2seed[kegg2seed['kegg']==kegg_id]['seed'].values[0]}
                    except IndexError: # there is no seed match for this kegg
                        mets_annotation[met] = { 'kegg.compound' : kegg_id}
        elif MET_INCH_COL in met_data.columns: # Specifically for INChIKey, important for different matchings
            if not pd.isna(met_data[MET_INCH_COL].values[0]):
                mets_annotation[met] = {'inchikey':met_data[MET_INCH_COL].values[0]}
        # WARNING: other annotations are not implemented
        else: # if no annotation, then fill it with pytfa like fake annotations
           mets_annotation[met] = {'seed_id' : 'fakeID_{}'.format(numerator)}
           numerator = numerator +1
           
        if MET_SMI_COL in met_data.columns: # For thermodynamic data search
            if not pd.isna(met_data[MET_INCH_COL].values[0]):
                mets_annotation[met].update({'smiles':met_data[MET_SMI_COL].values[0]})
        
                    
    return mets_formula, mets_name, mets_charge, mets_annotation

    
def merge_pthw_info(pthw_list):
    '''
    Merging the information of all the pathways to produce a target compound into
    a single pathway (hyper-pathway).

    '''
    pthw_id = list()
    rxns = list() 
    boundary_mets = list()
    bchd_pthws = list()
    
    for _, pthw in pthw_list.iterrows():
        this_rxns, target, this_boundary_mets, this_pthw_id, this_bchd_pthws = pthw_parser(pthw)
        rxns = rxns + [x for x in this_rxns if x not in rxns] 
        boundary_mets = boundary_mets + [x for x in this_boundary_mets \
                                                   if x not in boundary_mets]
        bchd_pthws = [x for x in this_bchd_pthws if x not in bchd_pthws]
        pthw_id = pthw_id + this_pthw_id
    
    target = target
    return rxns, target, boundary_mets, pthw_id, bchd_pthws

def safe_split(string, rep_1, rep_2='None', rep_3='None'):
    '''
    It splits string with rep_1 if possible, else with rep_2  

    Parameters
    ----------
    string : str
    rep_1 : str
    rep_2 : str

    Returns
    -------

    '''
    
    if rep_1 in string:
        return string.split(rep_1)
    elif rep_2 in string:
        return string.split(rep_2)
    elif rep_3 in string:
        return string.split(rep_3)
    else:
        raise Exception('Invalid format! Please, use either {} or {} or {} to split entries in a field').format(
            rep_1, rep_2, rep_3)

def find_pthw(pthw_ids, pthw_list):
    '''
    Based on the pthw_ids find the pathways that should be added among the pthw_list

    Parameters
    ----------
    bchd_pthws : list
    bchd_pthw_list : pd.DataFrame

    Returns
    -------
    pd.DataFrame

    '''

    pthw_dict = dict()
    
    for pthw_id in pthw_ids:
        pthw_dict.update(
            pthw_list[pthw_list['P_PR_UID']==pthw_id].to_dict(orient='index')
            )
            
            
    return pd.DataFrame.from_dict(pthw_dict, orient='index')


def import_bch_data(the_path):
    bchd_pthw_list = pd.read_csv(pjoin(the_path,'pathways_branching.tsv'), sep='\t')
    return bchd_pthw_list