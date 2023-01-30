__author__ = 'anastasia'

import os
import pandas as pd

def count_match(row, enzymes, df_bridgit):
    rxn_uids = [int(i) for i in row['Unnamed: 0'].replace('[','').replace(']','').replace("'",'').replace(" ",'').split(',')]
    ecs = df_bridgit[df_bridgit['rxnUID'].isin(rxn_uids)]['EC'].to_list()
    row['Num 4l EC matched'] = len(set(ecs).intersection(set(enzymes)))
    row['4l EC matched'] = ','.join(set(ecs).intersection(set(enzymes)))
    row['4l EC not matched'] = ','.join(set(enzymes) - set(ecs))
    return row

def count_3l_match(row, enzymes, df_bridgit):
    rxn_uids = [int(i) for i in row['Unnamed: 0'].replace('[','').replace(']','').replace("'",'').replace(" ",'').split(',')]
    ecs = ['.'.join(i.split('.')[0:3]) for i in df_bridgit[df_bridgit['rxnUID'].isin(rxn_uids)]['EC'].to_list()]
    row['Num 3l EC matched'] = len(set(ecs).intersection(set(enzymes)))
    row['3l EC matched'] = ','.join(set(ecs).intersection(set(enzymes)))
    row['3l EC not matched'] = ','.join(set(enzymes) - set(ecs))
    return row

def annotate_with_bridgit_match():
    df_bridgit = pd.read_csv(bridgit_results)
    for compound in compounds:
        print(compound)
        with open(folder_native_enzymes+compound) as f:
            enzymes = [i.strip()for i in f.readlines()]
        three_levels = ['.'.join(i.split('.')[0:3]) for i in enzymes] # get 3 level ecs
        print('All enzymes:', set(enzymes))
        print('All enzymes:', len(set(enzymes)))
        print('Three level enzymes:', set(three_levels))
        print('Three level enzymes:', len(set(three_levels)))
        df = pd.read_csv(folder_results_optimization+compound+'/final_scores.csv')
        df = df.apply(count_match, args = (enzymes, df_bridgit, ), axis=1)
        df = df.apply(count_3l_match, args = (three_levels, df_bridgit,), axis=1)
        df.to_csv(folder_results_optimization+compound+'/EC_match_natural.csv', index=False)

def get_overview_per_compound():
    for compound in compounds:
        print(compound)
        with open(folder_native_enzymes+compound) as f:
            enzymes = [i.strip()for i in f.readlines()]
        df_comp = pd.read_csv(folder_results_optimization+compound+'/EC_match_natural.csv')
        print('Min PW size', min(df_comp['Size']), 'Max PW size', max(df_comp['Size']))
        print('Min molar yield', min(df_comp['Molar Yield']), 'Max molar yield', max(df_comp['Molar Yield']))
        print('Min Num 4l EC matched', min(df_comp['Num 4l EC matched']), 'Max Num 4l EC matched', max(df_comp['Num 4l EC matched']))
        print('Min Num 3l EC matched', min(df_comp['Num 3l EC matched']), 'Max Num 3l EC matched', max(df_comp['Num 3l EC matched']))
        dict_match = dict()
        for enzyme in enzymes:
            dict_match[enzyme] = count_matches_enzyme(enzyme, df_comp)
        with open(folder_results_optimization+compound+'/enzyme_match_count.csv', 'w') as w:
            for enz, count in dict_match.items():
                w.write('{},{}\n'.format(enz, count))


def count_matches_enzyme(enzyme, df_comp):
    count = 0
    for index, row in df_comp.iterrows():
        if type(row['4l EC matched']) != float:
            if enzyme in row['4l EC matched'].split(','):
                count +=1
    return count

folder_results_optimization = "data/results_optmization/"
folder_native_enzymes = "data/8_compounds_native_enzymes/"
bridgit_results = '../data/ARBRE/bridgit_predictions.csv'

compounds = [i for i in os.listdir(folder_results_optimization) if not i.startswith('.')] # exclude .DS_store

annotate_with_bridgit_match()
get_overview_per_compound()



