#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 20:25:53 2021

@author: omid
"""
import pandas as pd
from os.path import join as pjoin


kegg2seed = pd.read_excel('../data/KEGG2SEED_update.xlsx',
                             sheet_name='Sheet1',
                             header=0)

def input_parser_netw(the_path):
    # a function to parse the pathways, metabolites and reactions input files
    # target: str, the name of taget compound
    
    met_list = pd.read_csv(pjoin(the_path,'compounds.tsv'), sep='\t')
    rxn_list = pd.read_csv(pjoin(the_path,'reactions.tsv'), sep='\t')
    # the metabolite IDs are numbers, we should convert them to strings
    met_list['M_PR_UID'] = met_list['M_PR_UID'].astype(str)
    rxn_list['R_PR_UID'] = rxn_list['R_PR_UID'].astype(str)
    
    return met_list, rxn_list

def input_parser_pthw(the_path):
    # a function to parse the pathways, metabolites and reactions input files
    # target: str, the name of taget compound
    
    met_list = pd.read_csv(pjoin(the_path,'compounds.tsv'), sep='\t')
    rxn_list = pd.read_csv(pjoin(the_path,'reactions.tsv'), sep='\t')
    pthw_list = pd.read_csv(pjoin(the_path,'pathways.tsv'), sep='\t')
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

def rxn_parser(rxns, rxn_list):
    mets = [] # list
    rxns_stoich = dict() # dict of dict
    remove_candid = []
    for rxn in rxns:
        rxn_dict = dict()
        rxn_data = rxn_list[rxn_list['R_PR_UID']==rxn]
        rxn_print = rxn_data['R_PR_STOICH'].values[0]
        participants = safe_split(rxn_print, ' ==> ', ' <=> ') # splitting reactants and products
        reactants = participants[0].split(' + ') # splitting reactants
        products = participants[1].split(' + ') # splitting products
        reactant_dict = {str(couple.split(' ')[1]):-1*float(couple.split(' ')[0])\
                         for couple in reactants}
        product_dict = {str(couple.split(' ')[1]):1*float(couple.split(' ')[0])\
                        for couple in products}
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
    mets_formula = dict() # dict
    mets_name = dict() # dict
    mets_charge = dict() # dict
    mets_annotation = dict() # dict
    numerator = 0 # to enumerate metabolites without annotations
    for met in mets:
        met_data = met_list[met_list['M_PR_UID']==met]
        mets_formula[met] = met_data['M_PR_FORMULA'].values[0] \
            if not pd.isna(met_data['M_PR_FORMULA'].values[0]) else ''
        mets_name[met] = met_data['M_PR_NAME'].values[0] if not pd.isna(met_data['M_PR_NAME'].values[0]) \
            else 'generic_metabolite'
        mets_charge[met] = met_data['M_PR_CHARGE'].values[0] if not pd.isna(met_data['M_PR_CHARGE'].values[0]) \
            else 0
        kegg_id = met_data['M_XR_KEGG'].values[0]
        other_fields = met_data['M_XR_ALL_OTHER'].values[0]
        if host == 'ecoli':
            if not pd.isna(kegg_id):
                try:
                    mets_annotation[met] = \
                        kegg2seed[kegg2seed['kegg']==kegg_id]['seed'].values[0]
                except IndexError: # there is no seed match for this kegg
                    mets_annotation[met] = ''
                    numerator = numerator +1
            else:
                mets_annotation[met] = ''
                numerator = numerator +1
        elif host == 'yeast':
            if not pd.isna(kegg_id):
                if kegg_id == 'C01328' or kegg_id == 'C00080': # for H2O and H+, use seed
                    mets_annotation[met] = \
                        kegg2seed[kegg2seed['kegg']==kegg_id]['seed'].values[0]
                else:
                    mets_annotation[met] = kegg_id
            else:
                # kegg id not found, look for chebi
                mets_annotation[met] = extract_chebi(other_fields, numerator)
                numerator = numerator +1
    return mets_formula, mets_name, mets_charge, mets_annotation

def extract_chebi(string, numerator):
    fields = string.split(',')
    chebi_ids = [field for field in fields if 'ChEBI' in field]
    if len(chebi_ids) == 0:
        annotation = ''
    else:
        annotation = chebi_ids[0].upper().replace(': ',':').split(' ')[0]
    return annotation

    
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

def safe_split(string, rep_1, rep_2):
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
    else:
        raise Exception('Invalid format! Please, use either {} or {} to split entries in a field').format(
            rep_1, rep_2)

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