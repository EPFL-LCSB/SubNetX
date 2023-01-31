#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 13:03:59 2021

@author: omid
"""
from connect_to_DB import closeConnection, fetchResults, connectToDatabase

from pytfa.io.json import load_json_model

import pandas as pd
from tqdm import tqdm


LCSBID = 'LCSBID'
NaN_ID = '<NA>'

model = load_json_model('iJO1366.json')

lexicon = pd.read_csv('LCSBids_iJO1366.csv', index_col = 0)
metabolite_bigg2lcsb = lexicon.astype('Int64').astype(str)
met_id_dict = metabolite_bigg2lcsb.to_dict(orient='index')

rxn_id_dict = dict()

for rxn in tqdm(model.reactions, desc='GEM reactions'):
    met_with_na_lcsb = [met.id for met in rxn.metabolites \
                        if met_id_dict[met.id][LCSBID] == NaN_ID]
    if len(met_with_na_lcsb) > 0: # only should search if all the metabolites have LCSB ID
        rxn_id_dict[rxn.id] = []
        continue
    
    if len(rxn.metabolites) == 1: # not for the exchange reactions
        rxn_id_dict[rxn.id] = []
        continue
    # We should prepare two inputs: the number of metabolites and their LCSB ID cancatenated in a string
    met_lsb_id = [met_id_dict[met.id][LCSBID] for met in rxn.metabolites]
    input_1 = ','.join(met_lsb_id)
    input_2 = str(len(rxn.metabolites))
    
    mydb = connectToDatabase()
    # match_rxn = fetchResults_comps_2_rxn(mydb,"1469424929,1467866080,1467866192,1467879801","4")
    match_rxn = fetchResults(mydb, input_1, input_2)
    mydb.commit()
    closeConnection(mydb)
    rxn_id_dict[rxn.id] = match_rxn
    print (match_rxn)
    
new = {}
for k,v in rxn_id_dict.items():
    new[k] = [str(x[0]) for x in v]
final_sol = {k:' | '.join(v) for k,v in new.items()}
data = pd.DataFrame(final_sol.items())
data.to_csv('rxn_LCSBID_iJO1366.csv')