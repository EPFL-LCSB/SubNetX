#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from pytfa.optim.utils import symbol_sum
from subnetx.optim.constraints import IntegerCutConstraint


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
    expr = symbol_sum(non_active_rxns) + symbol_sum(active_rxns)
    model.add_constraint(kind = IntegerCutConstraint, 
                            hook = model, 
                            expr = expr,
                            id_ = str(index),
                            lb = 1)


