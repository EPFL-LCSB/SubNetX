# -*- coding: utf-8 -*-
"""
Created on Thu May 29 18:29:35 2025

@author: oftad
"""
from pytfa.thermo.std import TEMPERATURE_0
from pytfa.thermo.equilibrator import compound_to_entry
import os
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_api import ComponentContribution, Q_

def append_thermodbs(thermodb1, thermodb2):
    '''
    It appends two thermo databases to retuen a single appended one.
    In case of conflicts, the second thermodb is favored. 

    Parameters
    ----------
    thermodb1 : a pytfa thermodb
    thermodb2 : a pytfa thermodb
    

    Returns
    -------
    big_thermodb : TYPE
        DESCRIPTION.

    '''
    
    if thermodb1['units'] != thermodb2['units']:
        raise ValueError("The units of the two databases must be the same.")
        
    big_thermodb = {}
    # the same name, units and cues as the second db
    big_thermodb["name"] = thermodb2["name"]
    big_thermodb["units"] = thermodb2["units"]
    big_thermodb["cues"] = thermodb2["cues"].copy()
    
    # start with the first db but overwrites with the second
    big_thermodb["metabolites"] = thermodb1["metabolites"].copy()
    big_thermodb["metabolites"].update(thermodb2["metabolites"])
    
    
    return big_thermodb

def build_thermo_novel(compound_df, met_smiles_dict, T=TEMPERATURE_0):
    """
    A function to generate thermo data for compounds not already included in the cache

    Parameters
    ----------
    compound_df : a pd.DataFrame
        DESCRIPTION.
    met_smiles_dict: dictionary with keys being metabolite IDs and values being SMILEs
    T : TYPE, optional
        DESCRIPTION. The default is TEMPERATURE_0.

    Returns
    -------
    thermo_data : TYPE
        DESCRIPTION.

    """
    lc = LocalCompoundCache()
    # from equilibrator_cache import create_compound_cache_from_zenodo
    # ccache = create_compound_cache_from_zenodo()
    if not os.path.exists('compounds.sqlite'):
        lc.generate_local_cache_from_default_zenodo('compounds.sqlite')
    lc.load_cache('compounds.sqlite')
    
    # no need to call add function as get_compound will add if new
    # lc.add_compounds(compound_df, mol_format="smiles") 
    cpd_results = {k:lc.get_compounds([v])[0].compound \
                   for k,v in met_smiles_dict.items()}

    # create pytfa-friendly thermodb
    cc = ComponentContribution()
    cc.temperature=Q_(str(T) + "K")

    thermo_data = {"name": "eQuilibrator", "units": "kJ/mol", "cues": {}}
    thermo_data["metabolites"] = {
        id_:compound_to_entry(id_,met, cc) for id_, met in cpd_results.items()
    }
    
    
    return thermo_data