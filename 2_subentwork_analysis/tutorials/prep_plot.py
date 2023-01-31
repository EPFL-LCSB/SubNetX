#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 18:11:10 2022

@author: omid
"""
import os

import pandas as pd
from pandas import ExcelWriter

from eval_netws import target_list, path_mod


def find_bar_labels(data_dict):
    # find the data labels that should be cobmined to plot 3-4-5 bars at the end
    
    if len(data_dict.keys()) < 5: # the number of bars is already small enough to plot
        label_dict = {str(round(k,2)):[k] for k in data_dict.keys()}
    else:
        tot_size = 0
        for val in data_dict.values():
            tot_size += len(val) # find the total number of pathways
        q_3 = round(tot_size/3)
        q_4 = round(tot_size/4)
        q_5 = round(tot_size/5)
        
        label_dict = {}
        
        this_label = 0
        labels = []
        for k,v in data_dict.items():
            this_label += len(v)
            labels += [k]
            if this_label > q_3: # the size of each group should not exceed 33% of the data, unless they all have the same key
                # the new key is big enough to be form one group on its own
                label_dict[str(round(k,2))] = [k]
                this_label -= len(v)
                labels.remove(k)
                if this_label != 0: # there was aome data before the latest key
                    new_key = '{}-{}'.format(round(min(labels),2),round(max(labels),2)) \
                        if round(min(labels),2)!=round(max(labels),2) else str(round(min(labels),2)) # use the average of grouped keys as a new key
                    label_dict[new_key] = labels
                    this_label = 0
                    labels = []
            elif this_label > q_5: # the size of each group is at least 20% of the data
                new_key = '{}-{}'.format(round(min(labels),2),round(max(labels),2)) \
                        if round(min(labels),2)!=round(max(labels),2) else str(round(min(labels),2)) # use the average of grouped keys as a new key
                label_dict[new_key] = labels
                this_label = 0
                labels = []
        
        if this_label != 0: # there is some data from the last group that is not yet grouped
            new_key = '{}-{}'.format(round(min(labels),2),round(max(labels),2)) \
                        if round(min(labels),2)!=round(max(labels),2) else str(round(min(labels),2)) # use the average of grouped keys as a new key
            label_dict[new_key] = labels
            this_label = 0
            labels = []         
        
    
    return label_dict


if __name__ == "__main__":
                     
    
    for target in target_list:
        ### preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        writer = ExcelWriter('{}/plot_tables.xlsx'.format(the_path)) 
        table = pd.read_csv('{}/{}/final_scores.csv'.format(path_mod,target))
        columns = ['Index','Size', 'Molar_Yield', 'Gram_Yield', 'BridgIT_score','BridgIT_weight', 'Thermodynamic_Eval', 'BridgIT_weight_2']
        table.columns = columns
        
        max_yield = table.Gram_Yield.max() # Universasl maximum yield
        min_size = table.Size.min() # Universasl minimum size
        
        ## group by size to find the highest yield, percent of thermodynamic feasibility, average BridgIT score
        size_dict = dict(tuple(table.groupby('Size')))
        result_dict = dict()
        
        # highest yield
        result_dict['size_yield'] = {}#{0:{'y':0}}
        for size, data in size_dict.items():
            result_dict['size_yield'][size] = {'y':data.Gram_Yield.max()}
            if data.Gram_Yield.max() == max_yield: # if we reach the maximum yield, there is no point to continue
                break
            
        # highest yield with feasibility check
        result_dict['size_yield_feas'] = {}#{0:{'y':0}}
        for size, d in size_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['size_yield_feas'][size] = {'y':data.Gram_Yield.max()}
            if data.Gram_Yield.max() == max_yield: # if we reach the maximum yield, there is no point to continue
                break
            
        # average BridgIT score
        result_dict['size_score'] = {}#{0:{'y':0}}
        for size, data in size_dict.items():
            result_dict['size_score'][size] = {'y':data.BridgIT_score.max()/size}
            
        # average BridgIT score with feasibility check
        result_dict['size_score_feas'] = {}#{0:{'y':0}}
        for size, d in size_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['size_score_feas'][size] = {'y':data.BridgIT_score.max()/size}
         
        # minimum BridgIT weight
        result_dict['size_weight'] = {}#{0:{'y':0}}
        for size, data in size_dict.items():
            result_dict['size_weight'][size] = {'y':data.BridgIT_weight.min()}
            
        # minimum BridgIT weight with feasibility check
        result_dict['size_weight_feas'] = {}#{0:{'y':0}}
        for size, d in size_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['size_weight_feas'][size] = {'y':data.BridgIT_weight.min()}
            
        # minimum BridgIT weight
        result_dict['size_weight_2'] = {}#{0:{'y':0}}
        for size, data in size_dict.items():
            result_dict['size_weight_2'][size] = {'y':data.BridgIT_weight_2.min()}
            
        # minimum BridgIT weight with feasibility check
        result_dict['size_weight_feas_2'] = {}#{0:{'y':0}}
        for size, d in size_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['size_weight_feas_2'][size] = {'y':data.BridgIT_weight_2.min()}
            
        # percent of feasibility
        result_dict['size_feasibility'] = {}
        # a preprocess to combine the keys that have small data
        label_dict = find_bar_labels(size_dict)
        for key, labels in label_dict.items():
            tot_size = 0
            tot_feas = 0
            for l in labels:
                data = size_dict[l]
                feas_dict = dict(tuple(data.groupby('Thermodynamic_Eval')))
                tot_size += len(data)
                try:
                    tot_feas += len(feas_dict['Feasible'])
                except KeyError:
                    pass
            result_dict['size_feasibility'][key] = {'y':tot_feas/tot_size}    
        
        
        
        ## group by yield to find the percent of infeasibility, average bridgit score, and minimum size
        yield_dict = dict(tuple(table.groupby('Molar_Yield')))
        
        
        # lowest size
        result_dict['yield_size'] = {} #{0:{'y':table.Size.min()}}
        for y, data in yield_dict.items():
            result_dict['yield_size'][y] = {'y':data.Size.min()}    
        # a postprocess to make sure size is increasing if yield increases 
        check = True
        while check:
            rem_key = []
            for index, (key,value) in enumerate(result_dict['yield_size'].items()):
                if index+1 != len(result_dict['yield_size']): # exclude the last index
                    if list(result_dict['yield_size'].values())[index]['y'] > list(result_dict['yield_size'].values())[index+1]['y']:
                        # it is ordered with increasing yield, so size should be also increasing
                        rem_key += [key]
            result_dict['yield_size'] = {k:v for k,v in result_dict['yield_size'].items() \
                                         if k not in rem_key}
            if len(rem_key) == 0:
                check = False
            
        # lowest size with feasibility check
        result_dict['yield_size_feas'] = {} #{0:{'y':table.Size.min()}}
        for y, d in yield_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['yield_size_feas'][y] = {'y':data.Size.min()}
        # a postprocess to make sure size is increasing if yield increases 
        check = True
        while check:
            rem_key = []
            for index, (key,value) in enumerate(result_dict['yield_size_feas'].items()):
                if index+1 != len(result_dict['yield_size_feas']): # exclude the last index
                    if list(result_dict['yield_size_feas'].values())[index]['y'] > list(result_dict['yield_size_feas'].values())[index+1]['y']:
                        # it is ordered with increasing yield, so size should be also increasing
                        rem_key += [key]
            result_dict['yield_size_feas'] = {k:v for k,v in result_dict['yield_size_feas'].items() \
                                         if k not in rem_key}
            if len(rem_key) == 0:
                check = False
        
        
        # average BridgIT score
        result_dict['yield_score'] = {} #{0:{'y':0}}
        for y, data in yield_dict.items():
            result_dict['yield_score'][y] = {'y':data.BridgIT_score.div(data.Size).min()}
            
        # average BridgIT score with feasibility check
        result_dict['yield_score_feas'] = {} #{0:{'y':0}}
        for y, d in yield_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['yield_score_feas'][y] = {'y':data.BridgIT_score.div(data.Size).min()}
        
        # minimum BridgIT weight
        result_dict['yield_weight'] = {} #{0:{'y':0}}
        for y, data in yield_dict.items():
            result_dict['yield_weight'][y] = {'y':data.BridgIT_weight.min()}
            
        # minimum BridgIT weight with feasibility check
        result_dict['yield_weight_feas'] = {} #{0:{'y':0}}
        for y, d in yield_dict.items():
            try:
                data = dict(tuple(d.groupby('Thermodynamic_Eval')))['Feasible']
            except KeyError:
                continue
            result_dict['yield_weight_feas'][y] = {'y':data.BridgIT_weight.min()}
            
        # percent of feasibility
        result_dict['yield_feasibility'] = {}
        # a preprocess to combine the keys that have small data
        label_dict = find_bar_labels(yield_dict)
        for key, labels in label_dict.items():
            tot_size = 0
            tot_feas = 0
            for l in labels:
                data = yield_dict[l]
                feas_dict = dict(tuple(data.groupby('Thermodynamic_Eval')))
                tot_size += len(data)
                try:
                    tot_feas += len(feas_dict['Feasible'])
                except KeyError:
                    pass
            result_dict['yield_feasibility'][key] = {'y':tot_feas/tot_size}    
                
                
        for key,data in result_dict.items():
            df = pd.DataFrame.from_dict(data, orient = 'index')
            df.to_excel(writer, sheet_name=key)
        writer.save()
            
            