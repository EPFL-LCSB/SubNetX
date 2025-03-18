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
                # 'ajmalicine', 
                # 'benzyl_cinnamate', 
                #    'benzylbenzoate', 
                #    'berberine', 
                # 'N_cinnamoyl_serotonin', 
                # 'Quercetin_3_O_6_acetylglucoside',
                # 'scopolamine', 
                # 'strictosidine',
                # 'tadalafil',
                # 'Biotin',
                # 'Ascorbate',
                # 'Azomycin',
                # 'Benzoate',
                # 'Caffeine',
                # 'Daidzein',
                # 'Kaempferitrin',
                # 'Aminopyrrolnitrin',
                # 'Genistein',
                # 'Nerolidol',
                # 'salicylamide',
                # '(+)-T-Muurolol',
                # "3,3'-diindolylmethane",
                # '4-Hydroxyphenylethanol',
                # '9-deacetoxyfumigaclavine_C',
                # '(-)-cyclopenine',
                # '(-)-Menthol',
                # 'alpha,alpha-Trehalose',
                # 'Nicotianamine',
                # 'beta-Nitropropanoate',
                # 'delta-Tocotrienol',
                # 'Cytidine',
                # 'Methyl_(indol-3-yl)acetate',
                # 'tetramethylpyrazine',
                # 'toyocamycin',
                # 'dehydrocyclopeptine',
                # 'Cephalosporin_C',
                # 'beta-Elemene',
                # 'Isopentenyl_adenosine',
                # 'Quinolinate',
                # 'raspberry_ketone',
                # 'Picropodophyllin',
                # 'DL-Cycloserine',
                # 'Pyrrolnitrin',
                # 'desferrioxamine_E',
                # 'curcumin',
                # 'hydroxytyrosol',
                # 'Resveratrol',
                # 'Violacein',
                # 'Artemisinate',
                # 'beta-Carotene',
                # 'Naringenin',
                # 'papaverine',
                # 'ferulic_acid',
                # 'Catechol',
                # 'E-cinnamate',
                # 'Gallic_acid',
                # '4-coumarate',
                # 'norcoclaurine',
                # 'dopamine',
                # 'chavicol',
                # 'Cynaroside',
                # 'Gastrodin',
                # 'reticuline',
                # 'rosmarinic_acid',
                # 'melatonin',
                # 'Mandelic_acid',
                # 'homoeriodictyol',
                # 'cis,cis-muconic_acid',
                # 'myricetin',
                # 'Sakuranetin',
                '2-phenylethylamine',
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