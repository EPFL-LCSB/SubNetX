#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:31:06 2020

@author: omid
"""


from add_subnet import model_builder

import os.path

organism = 'yeast'


path_mod = 'results_processed'

# target_list = pd.read_excel('../data/Compound-precursor.xlsx',
#                              sheet_name='Sheet1',
#                              header=0)['Compound'].tolist()
target_list = [
                'Aklavinone'
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
        model = model_builder(the_path, target, organism)