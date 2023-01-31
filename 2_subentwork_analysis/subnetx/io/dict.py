# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: Expression and thermodynamics-based models

.. moduleauthor:: Pierre Salvy

Make the model serializable
"""
from collections import OrderedDict, defaultdict

from tqdm import tqdm

import cobra.io.dict as cbd
from cobra.exceptions import SolverNotFound
from optlang.util import expr_to_json, parse_expr
from pytfa.io.dict import get_solver_string, var_to_dict, cons_to_dict, \
    obj_to_dict, rebuild_obj_from_dict
from pytfa.thermo.tmodel import ThermoModel
from pytfa.optim.variables import ReactionVariable, ModelVariable
from pytfa.optim.constraints import ReactionConstraint, ModelConstraint
from pytfa.optim.utils import get_all_subclasses

from ..core.chassis import ChassisModel
from ..optim import constraints, variables


SOLVER_DICT = {
    'optlang.gurobi_interface':'optlang-gurobi',
    'optlang.cplex_interface':'optlang-cplex',
    'optlang.glpk_interface':'optlang-glpk',
}

def make_subclasses_dict(cls):
    the_dict = {x.__name__:x for x in get_all_subclasses(cls)}
    the_dict[cls.__name__] = cls
    return the_dict

REACTION_VARIABLE_SUBCLASSES    = make_subclasses_dict(ReactionVariable)
REACTION_CONSTRAINT_SUBCLASSES  = make_subclasses_dict(ReactionConstraint)
MODEL_VARIABLE_SUBCLASSES       = make_subclasses_dict(ModelVariable)
MODEL_CONSTRAINT_SUBCLASSES     = make_subclasses_dict(ModelConstraint)

def metabolite_thermo_to_dict(metthermo):
    return metthermo.thermo.__dict__


def archive_variables(var_dict):
    obj = OrderedDict()

    obj['variables'] = []
    for classname,var in var_dict.items():
        obj[classname] = list(map(var_to_dict, var))

    return obj

def archive_constraints(cons_dict):
    obj = OrderedDict()

    for classname,cons in cons_dict.items():
        obj[classname] = list(map(cons_to_dict, cons))

    return obj
    

def model_to_dict(model):
    """

    :param model:
    :return:
    """

    # Take advantage of cobra's dict serialization for metabolites and
    # reactions
    obj = cbd.model_to_dict(model)

    obj['solver'] = get_solver_string(model)
    obj['objective'] = obj_to_dict(model)

    # Copy variables, constraints
    # obj['var_dict'] = archive_variables(model._var_kinds)
    # obj['cons_dict'] = archive_constraints(model._cons_kinds)
    obj['variables'] = list(map(var_to_dict, model._var_dict.values()))
    obj['constraints'] = list(map(cons_to_dict, model._cons_dict.values()))

    is_thermo = False

    if isinstance(model, ThermoModel):
        obj['kind'] = 'ThermoModel'
        obj['thermo_data'] = model.thermo_data #it's a dict
        obj['name'] = model.name
        obj['temperature'] = model.TEMPERATURE
        obj['min_ph'] = model.MIN_pH
        obj['max_ph'] = model.MAX_pH
        is_thermo = True

        # Relaxation info
        try:
            obj['relaxation'] = model.relaxation
        except AttributeError:
            pass

    if isinstance(model, ChassisModel):

        # Convenience attributes
        obj['rxn_lexicon'] = model.rxn_lexicon
        obj['met_lexicon'] = model.met_lexicon
        obj['organism'] = model.organism
        obj['products'] = [x.id for x in model.products]
        
        obj['hetero_mets'] = model.hetero_mets
        obj['hetero_rxns'] = model.hetero_rxns

        obj['kind'] = 'ChassisModel'
        
        
    # Metabolite and Reaction-level cleanup
    for rxn_dict in obj['reactions']:
        rxn = model.reactions.get_by_id(rxn_dict['id'])

        if is_thermo:
            _add_thermo_reaction_info(rxn, rxn_dict)
            
    # Peptides and Thermo
    for met_dict in obj['metabolites']:
        the_met_id = met_dict['id']
        is_peptide = False

        if is_thermo and not is_peptide: # peptides have no thermo
            the_met = model.metabolites.get_by_id(the_met_id)
            _add_thermo_metabolite_info(the_met, rxn_dict)
            met_dict['kind'] = 'Metabolite'

    return obj

def _add_thermo_reaction_info(rxn, rxn_dict):
    if hasattr(rxn, 'thermo'):
        rxn_dict['thermo'] = rxn.thermo

def _add_thermo_metabolite_info(met, met_dict):
    if hasattr(met, 'thermo'):
        met_dict['thermo'] = metabolite_thermo_to_dict(met)



def model_from_dict(obj, solver=None):
    # Take advantage of cobra's serialization of mets and reactions
    new = cbd.model_from_dict(obj)

    if solver is not None:
        try:
            new.solver = solver
        except SolverNotFound as snf:
            raise snf
    else:
        try:
            new.solver = obj['solver']
        except KeyError:
            pass

    # Populate variables and constraints
    if obj['kind'] == 'ChassisModel':
        new = ThermoModel(thermo_data=obj['thermo_data'],
                          model=new,
                          name=obj['name'],
                          temperature=obj['temperature'],
                          min_ph=obj['min_ph'],
                          max_ph=obj['max_ph'])
        new = init_thermo_model_from_dict(new, obj)
        new = ChassisModel(new, inplace=True)
        new = init_chassis_model_from_dict(new, obj)
        

    new._push_queue()

    # Force update GPR info
    for rxn in new.reactions:
        rxn.gene_reaction_rule = rxn.gene_reaction_rule

    for the_var_dict in tqdm(obj['variables'], desc='rebuilding variables'):
        this_id = the_var_dict['id']
        classname = the_var_dict['kind']
        lb = the_var_dict['lb']
        ub = the_var_dict['ub']
        try: #Backward compat
            scaling_factor = the_var_dict['scaling_factor']
        except KeyError:
            scaling_factor = 1

        if classname in REACTION_VARIABLE_SUBCLASSES:
            hook = new.reactions.get_by_id(this_id)
            this_class = REACTION_VARIABLE_SUBCLASSES[classname]
            nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    ub=ub,
                                    lb=lb,
                                    scaling_factor=scaling_factor,
                                    queue=True)
        elif classname in MODEL_VARIABLE_SUBCLASSES:
            hook = new
            this_class = MODEL_VARIABLE_SUBCLASSES[classname]
            nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    id_=this_id,
                                    ub=ub,
                                    lb=lb,
                                    scaling_factor=scaling_factor,
                                    queue=True)
        else:
            raise TypeError(
                'Class {} serialization not handled yet' \
                    .format(classname))

    new._push_queue()

    variable_parse_dict = {x.name:x for x in new.variables}

    for the_cons_dict in tqdm(obj['constraints'], desc='rebuilding constraints'):
        this_id = the_cons_dict['id']
        classname = the_cons_dict['kind']
        new_expr = parse_expr(the_cons_dict['expression'],
                              local_dict = variable_parse_dict)
        # str_expr = the_cons_dict['expression']
        #
        # # Sympify the expression so that we can substitute variables afterwards
        # sym_expr = sympify(str_expr)
        #
        # subs_dict = {x:new.variables.get(x.name) for x in sym_expr.free_symbols}
        #
        # new_expr = sym_expr.subs(subs_dict)
        lb = the_cons_dict['lb']
        ub = the_cons_dict['ub']

        if classname in REACTION_CONSTRAINT_SUBCLASSES:
            hook = new.reactions.get_by_id(this_id)
            this_class = REACTION_CONSTRAINT_SUBCLASSES[classname]
            nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr,
                                      ub=ub,
                                      lb=lb,
                                      queue=True)
        elif classname in MODEL_CONSTRAINT_SUBCLASSES:
            hook = new
            this_class = MODEL_CONSTRAINT_SUBCLASSES[classname]
            nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr, id_=this_id,
                                      ub=ub,
                                      lb=lb,
                                      queue=True)
        else:
            raise TypeError('Class {} serialization not handled yet' \
                            .format(classname))

    try:
        rebuild_obj_from_dict(new, obj['objective'])
    except KeyError:
        pass

    new.repair()
    return new

def init_thermo_model_from_dict(new, obj):
    for rxn_dict in obj['reactions']:
        rxn = new.reactions.get_by_id(rxn_dict['id'])

        if 'thermo' in rxn_dict:
            rxn.thermo = rxn_dict['thermo']

    for met_dict in obj['metabolites']:
        met = new.metabolites.get_by_id(met_dict['id'])

        if 'thermo' in met_dict:
            new._prepare_metabolite(met)

    # Relaxation info
    try:
        new.relaxation = obj['relaxation']
    except KeyError:
        pass

    return new

def init_chassis_model_from_dict(new, obj):
    new.rxn_lexicon = obj['rxn_lexicon']
    new.met_lexicon = obj['met_lexicon'] 
    new.organism = obj['organism']
    new.products = [new.metabolites.get_by_id(x) for x in obj['products']]
    
    new.hetero_mets = obj['hetero_mets']
    new.hetero_rxns = obj['hetero_rxns']
    return new
