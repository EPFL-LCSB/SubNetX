# -*- coding: utf-8 -*-
"""
.. module:: SubNetX
   :platform: Unix, Windows

.. moduleauthor:: SubNetX team

Constraints declarations

"""

from pytfa.optim import ReactionConstraint, ModelConstraint


class UpperBoundCoupling(ReactionConstraint):
    """
    Class to represent a forward directionality coupling with binary variable
    Looks like :
    RU_rxn: F_rxn - R_rxn - M BFUSE_rxn < 0
    """

    prefix = 'RU_'
    
class LowerBoundCoupling(ReactionConstraint):
    """
    Class to represent a forward directionality coupling with binary variable
    Looks like :
    RL_rxn: F_rxn - R_rxn + M BFUSE_rxn > 0
    """

    prefix = 'RL_'
    
class IntegerCutConstraint(ModelConstraint):
    ''' A set of constraints to make the previously found solutions infeasible.
    '''
    prefix = 'ICC_'
    
class NonZeroForce(ReactionConstraint):
    """
    Class to to enforce that the absolute flux must be nonzero if the binary
    variable is 1:
    NZF_rxn: F_rxn + R_rxn - eps * BFUSE_rxn > 0
    """

    prefix = 'NZF_'