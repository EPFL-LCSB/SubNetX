#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:16:33 2021

@author: omid
"""
import pandas as pd
import os
from os.path import join as pjoin


path_mod = '../../results_processed-master/EPFL'

target_list = pd.read_excel('../data/Compound-precursor.xlsx',
                             sheet_name='Sheet1',
                             header=0)['Compound'].tolist()

def input_parser(the_path):
    # a function to parse the pathways, metabolites and reactions input files
    # target: str, the name of taget compound
    
    pthw_list = pd.read_csv(pjoin(the_path,'pathways.tsv'), sep='\t')
    target_col = pthw_list['P_PR_TARGET']
    target_id = str(target_col[0])
    return target_id
    
    

if __name__ == "__main__":
    ret = {}
    for target in target_list:
        # preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        ret[target] = input_parser(the_path)
        