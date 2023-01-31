#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 18:17:50 2022

@author: omid
"""
import os
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import math

from eval_netws import target_list, path_mod

# plt.style.use('seaborn-white')
# plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'

colorset = [
    'red',
    'black',
    'green',
    'blue',
    ]

def find_data_range(super_data):
    
    min_x = 1000
    max_x = 0
    min_y = 1000
    max_y = 0
    for _, data in super_data.items():
        data.columns = ['x', 'y']
        x = data['x'].tolist()
        y = data['y'].tolist()
        if min(x) < min_x:
            min_x = min(x)
        if max(x) > max_x:
            max_x = max(x)
        if min(y) < min_y:
            min_y = min(y)
        if max(y) > max_y:
            max_y = max(y)
    
    return min_x, max_x, min_y, max_y


def create_2d_plot(target, super_data, xlabel, ylabel, image_name,
                   xtype='float', ytype='int',
                   **kwargs):
    
    hfont = {'fontname':'Helvetica', 'fontsize':18}
    
    fig, axs = plt.subplots(1,len(super_data))
    fig.suptitle(target.capitalize().replace('_',' '),**hfont)
    
    # if it is only one figure the axs object is not indexable, fix this by creating a list
    try:
        _ = axs[0]
    except TypeError:
        axs = [axs]
        
    
    for ind, (key, data) in enumerate(super_data.items()):
        data.columns = ['x', 'y']
        # plt.plot(data['x'].tolist(),data['y'].tolist())
        
        x = data['x'].tolist()
        y = data['y'].tolist()
        axs[ind].plot(x,y, color=colorset[ind], 
                      **kwargs
                      )
    
        axs[ind].set_title(key, 
                           # y=-0.25, 
                           fontdict={'fontsize': 9, 'fontweight': 'medium'})
        
        # find the x and y labels
        # if more than one plot, try to standardize the axes based on highest variation
        min_x, max_x, min_y, max_y = find_data_range(super_data)
            
        if xtype == 'int':
            min_x_l = min_x
            max_x_l = max_x
            xsticks = np.arange(min_x_l, max_x_l+1, 1)
        elif xtype == 'float':
            # diff = max(x) - min(x)
            # n_d = [x for x in range(0,20) \
            #        if 10**x*diff==int(10**x*diff)][0]
            min_x_l = math.floor(min_x*10)/10 # one digit after decimal
            max_x_l = math.ceil(max_x*10)/10
            xsticks = np.linspace(min_x_l, max_x_l, 6, endpoint=True)
        axs[ind].set_xticks(xsticks)
        if ytype == 'int':
            min_y_l = min_y 
            max_y_l = max_y 
            ysticks = np.arange(min_y_l, max_y_l+1, 1)
        elif ytype == 'float':
            min_y_l = math.floor(min_y*10)/10 # one digit after decimal
            max_y_l = math.ceil(max_y*10)/10
            ysticks = np.linspace(min_y_l, max_y_l, 6, endpoint=True)
        axs[ind].set_yticks(ysticks)
        # the intersection
        # ax.spines['left'].set_position(('data', min_y_l))
        # ax.spines['bottom'].set_position(('data', min_x_l))
        # ax.spines['left'].set_bounds(min_y_l,max_y_l)
        # ax.spines['bottom'].set_bounds(min_x_l,max_x_l)
        for spine in ['top', 'right']:
            axs[ind].spines[spine].set_visible(False)
        
        #  naming the axes
        axs[ind].set_xlabel(xlabel, fontdict={'fontsize': 12, 'fontweight': 'medium'})
        axs[ind].set_ylabel(ylabel, fontdict={'fontsize': 12, 'fontweight': 'medium'})
    
    fig.set_figwidth(12)
    fig.set_figheight(4)
    
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #     ax.label_outer()    
    # function to show the plot 
    plt.show()
    
    image_format = 'svg' # e.g .png, .svg, etc.    
    fig.savefig(image_name, format=image_format, dpi=1200)
    
    return  
    

def create_bar_plot(target, super_data, xlabel, ylabel, image_name=None,
                    w_dict = {'width_ratios': [1, 1]}):
    
    hfont = {'fontname':'Helvetica', 'fontsize':18}
    
    fig, axs = plt.subplots(1,len(super_data),gridspec_kw=w_dict)
    fig.suptitle(target.capitalize().replace('_',' '),**hfont)
    
    # if it is only one figure the axs object is not indexable, fix this by creating a list
    try:
        _ = axs[0]
    except TypeError:
        axs = [axs]
        
        
    for ind, (key, data) in enumerate(super_data.items()):
        data.columns = ['x', 'y']
        
        x = data['x'].tolist()
        y = data['y'].tolist()
        axs[ind].bar(x,y,color='orange'
                      )
    
        axs[ind].set_title(key, 
                            y=-0.15, 
                           fontdict={'fontsize': 12, 'fontweight': 'medium'})
        #  naming the axes
        # axs[ind].set_xlabel(xlabel[ind], fontdict={'fontsize': 12, 'fontweight': 'medium'})
        if ind == 0: # only for the outer plot
            axs[ind].set_ylabel(ylabel, fontdict={'fontsize': 12, 'fontweight': 'medium'})
        
        
    fig.set_figwidth(9)
    fig.set_figheight(6)
    
    plt.show()
    
    image_format = 'svg' # e.g .png, .svg, etc.    
    fig.savefig(image_name, format=image_format, dpi=1200)
    
    
    return

if __name__ == "__main__":
                     
    
    for target in target_list:
        ### preliminaries
        target = target.replace(' ', '_')
        the_path = '{}/{}'.format(path_mod,target)
        print('The current compound is {}.'.format(target))
        if not os.path.isdir(the_path):
            print('{} is not found.'.format(target))
            continue
        
        # plot yield against size
        data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_size')
        data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_size_feas')
        data = {
                'without thermodynamics':data_1,
                'with thermodynamics'   :data_2
                }
        create_2d_plot(target, data, image_name = '{}/SizeVsYield.svg'.format(the_path),
                       xlabel='Yield (mol mol$^{-1}$)', ylabel='Size (number of the reactions)',
                       linestyle = '--', marker = 'o',
                       fillstyle='full', markersize=3)
        
        
        # plot size against yield
        # data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_yield')
        # data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_yield_feas')
        # data = {
        #         'without thermodynamics':data_1,
        #         'with thermodynamics'   :data_2
        #         }
        # create_2d_plot(target, data, image_name = '{}/YieldVsSize.svg'.format(the_path),
        #                xlabel='Size', ylabel='Yield',
        #                linestyle = '--', marker = 'x',
        #                xtype = 'int', ytype = 'float')
        
        
        
        # data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_score')
        # data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_score_feas')
        # data = {
        #         'without thermodynamics':data_1,
        #         'with thermodynamics'   :data_2
        #         }
        # create_2d_plot(target, data, image_name = '{}/ScoreVsSize.svg'.format(the_path),
        #                xlabel='Size', ylabel='BridgIT Score',
        #                 linestyle = '--', marker = 'x',
        #                xtype = 'float', ytype = 'float')
        
        
        
        # data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_score')
        # data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_score_feas')
        # data = {
        #         'without thermodynamics':data_1,
        #         'with thermodynamics'   :data_2
        #         }
        # create_2d_plot(target, data, image_name = '{}/ScoreVsYield.svg'.format(the_path),
        #                xlabel='Yield', ylabel='BridgIT Score',
        #                 linestyle = '--', marker = 'x',
        #                xtype = 'float', ytype = 'float')
        
        
        
        # plot weight against yield
        # data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_weight')
        # data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_weight_feas')
        # data = {
        #         'without thermodynamics':data_1,
        #         'with thermodynamics'   :data_2
        #         }
        # create_2d_plot(target, data, image_name = '{}/WeightVsYield.svg'.format(the_path),
        #                xlabel='Yield', ylabel='BridgIT Weight',
        #                # linestyle = '--', marker = 'x',
        #                xtype = 'float', ytype = 'float')
        
        
        # plot weight against size
        data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_weight')
        data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_weight_feas')
        data = {
                'without thermodynamics':data_1,
                'with thermodynamics'   :data_2
                }
        create_2d_plot(target, data, image_name = '{}/WeightVsSize.svg'.format(the_path),
                       xlabel='Size', ylabel='BridgIT Weight',
                       linestyle = '--', marker = 'x',
                       xtype = 'int', ytype = 'float')
        
        # plot weight against size
        data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_weight_2')
        data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_weight_feas_2')
        data = {
                'without thermodynamics':data_1,
                'with thermodynamics'   :data_2
                }
        create_2d_plot(target, data, image_name = '{}/WeightVsSize_2.svg'.format(the_path),
                       xlabel='Size', ylabel='BridgIT Weight',
                       linestyle = '--', marker = 'x',
                       xtype = 'int', ytype = 'float')
        
        # plot feasibility ratio against size and yield
        data_1 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='size_feasibility')
        data_2 = pd.read_excel('{}/plot_tables.xlsx'.format(the_path), sheet_name='yield_feasibility')
        data = {
                'Size (number of the reactions)' :data_1,
                'Yield (mol mol$^{-1}$)':data_2
                }
        create_bar_plot(target, data, image_name = '{}/FeasibilityRatio.svg'.format(the_path),
                       xlabel='None', ylabel='Fraction of the feasible pathways',
                       w_dict = {'width_ratios': [1, 1.4]})
        
        