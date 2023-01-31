#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:31:06 2020

@author: omid
"""


from addSubnetEcoli import model_builder

import os.path



path_mod = '../tutorials/results_processed-master/EPFL'

# target_list = pd.read_excel('../data/Compound-precursor.xlsx',
#                              sheet_name='Sheet1',
#                              header=0)['Compound'].tolist()
target_list = [
                'ajmalicine', 
                'benzyl_cinnamate', 
                   'benzylbenzoate', 
                   'berberine', 
                'N_cinnamoyl_serotonin', 
                'Quercetin_3_O_6_acetylglucoside',
                'scopolamine', 
                'strictosidine',
                'tadalafil',
               ]

if __name__ == "__main__":
    for target in target_list:
        # preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        
        # Doing pathway evaluation :-)
        model = model_builder(the_path, target)