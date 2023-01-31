#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:31:06 2020

@author: omid
"""


from thermoDBconverter import matToPy
from thermoModelbuilder import model_builder

import pandas as pd
import os.path
from scipy.io import loadmat
import pickle
import zlib


path_mod = '../work/results_processed-master/EPFL'
thdb_mat = ['thermodata_ecoli.mat', 'thermodata_yeast.mat']
thdb_py = {'thermodata_ecoli.mat' : 'thermo_data_ecoli.thermodb',
           'thermodata_yeast.mat' : 'thermo_data_yeast.thermodb'}

target_list = pd.read_excel('../data/Compound-precursor.xlsx',
                             sheet_name='Sheet1',
                             header=0)['Compound'].tolist()

if __name__ == "__main__":
    for target in target_list:
        # preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        # converting thermodb.mat to thermodb.thermodb
        for org_thdb in thdb_mat:
            thdb_path = '{}/{}'.format(the_path,org_thdb)
            save_path = '{}/{}'.format(the_path,thdb_py[org_thdb])
            ReactionDB = loadmat(thdb_path)
            # Guess the key of the array to the data of the database
            print("Guessing the name of the database...")
            name = None
            for key in ReactionDB:
                if key[:2] == "__": # Name begins with "__" (private attribute)
                    continue
                if name == None:
                    name = key
                else:
                    print("Cannot guess the name of the database ! Found at least 2 candidates : ")
                    print(name + " and " + key)
                    break
        
            if name == None:
                print("Cannot find name, no candidates found")
                continue
        
            print("Found name : " + name)
            ReactionDB = ReactionDB[name]
        
            # Convert all the matlab data to a more python-friendly version
            (units, metabolites, cues) = matToPy(ReactionDB)
        
            # Finally, save the data
            try:
                with open(save_path, 'wb') as file:
                    # Convert it to string with Pickle
                    data = pickle.dumps({
                        "name":name,
                        "units": units,
                        "metabolites":metabolites,
                        "cues":cues
                    }, pickle.HIGHEST_PROTOCOL)
        
                    # Compress it with zlib
                    data = zlib.compress(data)
        
                    # And write it to disk
                    file.write(data)
            except:
                print('Could not save the data')
                continue
        
        # Doing pathway evaluation :-)
        model_builder(the_path, target, hosts = ['ecoli','yeast'])