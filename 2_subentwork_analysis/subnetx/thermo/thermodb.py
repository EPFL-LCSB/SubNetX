# -*- coding: utf-8 -*-
"""
Created on Thu May 29 18:29:35 2025

@author: oftad
"""


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