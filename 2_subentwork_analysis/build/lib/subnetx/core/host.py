#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:37:42 2020

@author: omid
"""
import pandas as pd

from cobra import Metabolite, Reaction

from cobra import Model
from pytfa.thermo.tmodel import ThermoModel

from ..io.parser import met_parser, rxn_parser


LCSBID = 'LCSBID'
ANNOTATION = 'seed_id'
EXPORTID = 'DM_'

def hash2id_met(mets, met_id_prefix): 
        # assigns a generic ID to the metabolites
        mets_id = {met : met_id_prefix+str(ind) for ind,met in enumerate(mets)}
        return mets_id
    
def hash2id_rxn(rxns, rxn_id_prefix):
    # assigns a generic ID to the reactions
    rxns_id = {rxn : rxn_id_prefix+str(ind) for ind,rxn in enumerate(rxns)}
    return rxns_id

class HostModel(ThermoModel):
    
    def __init__(self, model,
                 organism='', met_lexicon=None, rxn_lexicon=None, 
                 products=[], inplace=True,
                 *args, **kwargs):
        '''
        

        Parameters
        ----------
        model : cobra.model
        organism : str, optional
            The name of host organism, e.g. ecoli or yeast. The default is ''.
        lexicon : pd.DataFrame, optional
            The LCSB annotation for the metabolites that are present or added. The default is None.
        products : list, optional
            The list of target compounds that can be produced by this host. The default is [].



        Returns
        -------
        None.

        '''
        
        if not inplace:
            new = model.copy()
            self.__dict__ = new.__dict__
        else:
            # new = me_model
            self.__dict__ = model.__dict__

        self.organism = organism
        self._met_lexicon = met_lexicon
        self._rxn_lexicon = rxn_lexicon
        self._products = products
        
        self._hetero_mets = {} # a dict to keep mets_id
        self._hetero_rxns = [] # a list to keep rxns
        
    @property
    def met_lexicon(self):
        return self._met_lexicon
    
    @met_lexicon.setter
    def met_lexicon(self, value):
        self._met_lexicon = value
        
    @property
    def rxn_lexicon(self):
        return self._rxn_lexicon
    
    @rxn_lexicon.setter
    def rxn_lexicon(self, value):
        self._rxn_lexicon = value
        
    @property
    def products(self):
        return self._products
    
    @products.setter
    def products(self, value):
        self._products = value
        
    @property
    def hetero_mets(self):
        return self._hetero_mets
    
    @hetero_mets.setter
    def hetero_mets(self, value):
        self._hetero_mets = value
        
    @property
    def hetero_rxns(self):
        return self._hetero_rxns
    
    @hetero_rxns.setter
    def hetero_rxns(self, value):
        self._hetero_rxns = value
        
    def add_mets_rxns(self, met_list, rxn_list, 
                      met_id_prefix = 'metabolite_', rxn_id_prefix = 'reaction_'):
        # 1. To add reactions of the pathway or network
        rxns = self.hetero_rxns # The reactions for the pathway or network
        
        mets, rxns_stoich = rxn_parser(rxns, rxn_list)
        # self.hetero_mets = mets
        
        rxns_id = hash2id_rxn(rxns, rxn_id_prefix)
        #TODO: write a similar function for reactions
        # rxns_id = self._reconcile_rxns(rxns_id)
        
        reactions = [Reaction(rxns_id[id_],
                              lower_bound = -1000, 
                              upper_bound = 1000) for id_ in rxns \
                              if rxn_id_prefix in rxns_id[id_]] # only add new reactions
        
        self.add_reactions(reactions)
        
        # 2. To add metabolites of the pathway or network
        # mets = self.hetero_mets # The metabolites for the pathway or network
        
        mets_formula, mets_name, mets_charge, \
            mets_annotation = met_parser(mets, met_list, self.organism)
            
        mets_id = hash2id_met(mets, met_id_prefix)
        mets_id = self._reconcile_mets(mets_id)
        self.hetero_mets.update(mets_id) # The metabolites for the pathway or network

        metabolites = [Metabolite(mets_id[id_],
                                  formula=mets_formula[id_],
                                  name=mets_name[id_], 
                                  charge=int(mets_charge[id_]),
                                  compartment='c') for id_ in mets \
                                   if met_id_prefix in mets_id[id_]] # only add new metabolites
        
        self.add_metabolites(metabolites)
        for id_ in mets:
            if mets_id[id_] != '': # make sure to have a real annotation
                self.metabolites.get_by_id(mets_id[id_]).annotation = \
                    {ANNOTATION : mets_annotation[id_]} # this is used to search in thermodb
    
        ## This part connects metabolites and reactions
        # Updating the metabolites in the added reactions:
        for id_ in rxns:
            the_reaction = self.reactions.get_by_id(rxns_id[id_])
            the_reaction.add_metabolites({mets_id[k]:v for k,v in rxns_stoich[id_].items()})
            
        # Finally, update the lexicon based on the added metabolites, so it can be used when the function is called again
        self.update_met_lexicon()
        self.update_rxn_lexicon(rxns_id)
    
    def _reconcile_mets(self, mets_id):
    # associate pathway metabolites with GEM metabolites, if already exist
        lexicon = pd.DataFrame.from_dict(self.met_lexicon,  orient='index')
        for hash_ in mets_id.keys():
            try:
                mets_id[hash_] = [x for x in lexicon[lexicon[LCSBID]== \
                                              hash_].index
                                              if '_e' not in x][0] # any compartment other than extracellular is acceptable
            except IndexError:
                continue
        
        return mets_id
    
    def add_target(self, target):
        # To add the target of a pathway or network
        # target is the string for the meyabolite ID of the target
        mets_id = self.hetero_mets
        target_rxn = [Reaction(EXPORTID + mets_id[target],
                            lower_bound = 0, # only to be exported
                            upper_bound = 1000)] # export the target
        # it is possible that the target reaction is already in the model
        try:
            self.reactions.get_by_id(target_rxn[0].id)
        except KeyError:
            self.add_reactions(target_rxn)
            target_rxn[0].add_metabolites({mets_id[target]:-1})
            self.products += [self.metabolites.get_by_id(mets_id[target])]
        return target_rxn[0].id
        
    def add_boundary_rxns(self, boundary_mets):
        # To add exchange reactions for the boundary metabolites for a pathway
        mets_id = self.hetero_mets
        ex_rxns = [Reaction(EXPORTID + mets_id[id_],
                            lower_bound = 0, # only to be exported
                            upper_bound = 1000) for id_ in boundary_mets]
        
        self.add_reactions(ex_rxns)
        
        # Now, modify the stoichiometries and metabolite annotations
        
        for id_ in boundary_mets:
            self.reactions.get_by_id(EXPORTID + mets_id[id_]).add_metabolites({mets_id[id_]:-1})
            
    def update_met_lexicon(self):
    
        old_dict = self.met_lexicon
        mets_id  = self.hetero_mets
        new_dict =  {v:{LCSBID : k} for k,v in mets_id.items()}
        old_dict.update(new_dict)
        self.met_lexicon = old_dict
        return
        
    def update_rxn_lexicon(self, rxns_id):
    
        old_dict = self.rxn_lexicon
        new_dict =  {v:{LCSBID : k} for k,v in rxns_id.items()}
        old_dict.update(new_dict)
        self.rxn_lexicon = old_dict
        return