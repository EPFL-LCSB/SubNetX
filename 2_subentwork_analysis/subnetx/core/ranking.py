#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 11:00:56 2021

@author: omid
"""
from tqdm import tqdm

import numpy as np
from cobra import Reaction
from pytfa.optim.utils import symbol_sum
from pytfa.optim.variables import ForwardBackwardUseVariable
from subnetx.optim.constraints import IntegerCutConstraint, UpperBoundCoupling,\
    LowerBoundCoupling, NonZeroForce

bigM = 1000
epsilon = 0.01

EXPORTID = 'DM_'

def assign_bin_var2rxn(model, check_rxn_id = ''):
    
    for rxn in tqdm(model.reactions, desc='adding binary variables'):
        if check_rxn_id in rxn.id:
            U_rxn = model.add_variable(ForwardBackwardUseVariable, rxn)
            model.add_constraint(UpperBoundCoupling, rxn, rxn.flux_expression - bigM * U_rxn,
                                 ub=0)
            model.add_constraint(LowerBoundCoupling, rxn, rxn.flux_expression + bigM * U_rxn,
                                 lb=0) 
            
    model.repair()
    
def bias_passing_intermediate(model, intermediates, weight = None):
    # To bias the pathway to pass through certain metabolites: 
    if not hasattr(intermediates, '__iter__'):
        intermediates = [intermediates]
        
    if isinstance(intermediates[0], str):
        intermediates = [model.metabolites.get_by_id(x) for x in intermediates]
        
    b_use_vars = model.get_variables_of_type(ForwardBackwardUseVariable)
    model.objective = symbol_sum([var.variable for var in b_use_vars] )
    
    for this_intermediate in intermediates:
        for rxn in this_intermediate.reactions:
            # We should couple the binary variable to the flux so that if it is 1, the flux must be nonzero
            try:
                U_rxn = b_use_vars.get_by_id(rxn.id)
            except KeyError:
                continue
            model.add_constraint(NonZeroForce, rxn, 
                                 rxn.forward_variable +  rxn.reverse_variable\
                                     - epsilon * U_rxn,
                                 lb=0)
            # We should set a different penalty in the objective
            if weight is None:
                model.objective.set_linear_coefficients({U_rxn.variable:-1}) # these reactions are rewarded instead of being penalized
            elif isinstance(weight, dict):
                the_coef = weight[rxn.id]
                model.objective.set_linear_coefficients({U_rxn.variable:the_coef})
            elif isinstance(weight, int) or isinstance(weight, float):
                model.objective.set_linear_coefficients({U_rxn.variable:weight})
            
        # We should add a demand reaction to allow its pull-out
        ex_rxn = Reaction(EXPORTID + this_intermediate.id,
                    lower_bound = 0, # only to be exported
                    upper_bound = 1000)
        model.add_reactions([ex_rxn])
        ex_rxn.add_metabolites({this_intermediate:-1})
    
    
    model.objective_direction = 'min'
    return

def add_integer_cuts(model, varaibles, index=''):
    '''
    Adds integer cut constraints to enumerate diffrent solutions

    Parameters
    ----------
    model : TYPE
        DESCRIPTION.
    varaibles : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    active_rxns = []
    non_active_rxns = []
    for b_var in varaibles:
        try:
            var = b_var.variable
        except AttributeError: # it is an internal variable
            pass
        
        if np.isclose(var.primal, 0):
            non_active_rxns.append(var)
        elif np.isclose(var.primal, 1):
            active_rxns.append(1-var)
            
    # Adding the integer cut constraints to enumerate the other solutions
    expr = symbol_sum(active_rxns) # + symbol_sum(non_active_rxns) # if we do it only on active reactions, it's better
    model.add_constraint(kind = IntegerCutConstraint, 
                            hook = model, 
                            expr = expr,
                            id_ = str(index),
                            lb = 1)
