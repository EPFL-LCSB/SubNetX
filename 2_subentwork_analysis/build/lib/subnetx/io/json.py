# -*- coding: utf-8 -*-
"""
.. module:: subnetx
   :platform: Unix, Windows

.. moduleauthor:: Omid Oftadeh

JSON serialization
"""

import json

from .dict import model_from_dict, model_to_dict
from pytfa.io.json import MyEncoder, check_json_extension

def save_json_model(model, filepath):
    """
    Saves the model as a JSON file

    :param model:
    :param filepath:
    :return:
    """

    filepath = check_json_extension(filepath)
    obj = model_to_dict(model)

    with open(filepath, 'w') as fid:
        json.dump(obj, fid, cls=MyEncoder)


def load_json_model(filepath, solver = None):
    """
    Loads a model from a JSON file

    :param filepath:
    :param solver:
    :return:
    """

    filepath = check_json_extension(filepath)
    with open(filepath, 'r') as fid:
        obj = json.load(fid)

    model = model_from_dict(obj, solver=solver)
    return model