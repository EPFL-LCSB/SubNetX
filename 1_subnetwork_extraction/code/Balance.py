__author__ = 'anastasia'

import pandas as pd
from rdkit import Chem
import json
import os

class Balance():
    def __init__(self, data):
        self.data = data
        self.df_compounds = data.df_compounds
        self.df_reactions = data.df_reactions_comp
        self.get_compounds_atoms_dict() # load the dictionary with count of atom types per compound


    # creating atom composition dictionary for each compound
    def process_smiles(self, smiles):
        """
        This function processes smiles and gets dict with number of atoms per atom type
        :param smiles: input SMILES that have to be converted to the dict with number of atoms per atom type
        :return: dict with number of atoms per atom type
        """
        #print(smiles)
        molecule_atoms_dict = dict()
        if type(smiles) == float:
            return molecule_atoms_dict
        m = Chem.MolFromSmiles(smiles)
        if m is None: return molecule_atoms_dict
        for atom in m.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num not in molecule_atoms_dict.keys():
                molecule_atoms_dict[atomic_num]=0
            molecule_atoms_dict[atomic_num]+=1
        return molecule_atoms_dict

    def get_compounds_atoms_dict(self):
        """
        This function returns dict with
        :return:
        """
        reaction_network = self.data.data_folder+self.data.reaction_network

        if os.path.exists(reaction_network+'/compound_atoms.json') :
            with open(reaction_network+'/compound_atoms.json') as f:
                self.compound_atoms = json.load(f)
                #self.compound_atoms = self.convert_keys_to_int(self.compound_atoms)
            return

        # if the json with the atom composition of all the compounds does not exist yet, create it
        self.compound_atoms = dict() # creating dict to store the atom composition for every LCSB compound
        for index, row in self.df_compounds.iterrows():
            #print(index)
            smiles = row['SMILES']
            molecule_atoms_dict = self.process_smiles(smiles)
            if not molecule_atoms_dict:
                print(index)
                print(row['cUID'], row['COMMON_NAME'], row['SMILES'])
            self.compound_atoms[row['cUID']] = molecule_atoms_dict

        with open(reaction_network+'/compound_atoms.json', 'w') as f:
            json.dump(self.compound_atoms, f)

        return


    def convert_keys_to_int(self, d):
        """
        json dumps keys as strings, so they have to be converted back to int
        :param d: json dictionary
        :return:
        """
        new_dict = {}
        for k, v in d.items():
            try:
                new_key = int(k)
            except ValueError:
                new_key = k
            if type(v) == dict:
                v = self.convert_keys_to_int(v)
            new_dict[new_key] = v
        return new_dict

    #checking balance for reaction
    def expand_compounds_list_with_stoichimetry(self, comp_list_with_stoich):
        '''
        this funstion expands the list to include molecule as many times as it is mentioned in the stoich coef
        '''
        compounds_list = []
        for cmp_stoich in comp_list_with_stoich:
            stoich = int(cmp_stoich.split(' ')[0])
            comp = cmp_stoich.split(' ')[1]
            compounds_list.extend([comp]*abs(stoich))
        return compounds_list

    def check_balance(self, rxn):
        reactants_with_stoich = [i for i in rxn.split(';') if i.startswith('-')]
        products_with_stoich = [i for i in rxn.split(';') if not i.startswith('-')]

        reactants = self.expand_compounds_list_with_stoichimetry(reactants_with_stoich)
        products = self.expand_compounds_list_with_stoichimetry(products_with_stoich)

        total_atoms_reactants = dict()
        total_atoms_products = dict()

        for cmp in reactants:
            cmp_dict = self.compound_atoms[cmp]
            if not cmp_dict:
                return 'EXCLUDE INCLOMPLETE COMPOUND STRUCTURES REACTION'
            total_atoms_reactants = {k: total_atoms_reactants.get(k, 0) + cmp_dict.get(k, 0) for k in set(total_atoms_reactants) | set(cmp_dict)} #adding atoms one by one molecule to the total atoms dictionary

        for cmp in products:
            cmp_dict = self.compound_atoms[cmp]
            if not cmp_dict:
                return 'EXCLUDE INCLOMPLETE COMPOUND STRUCTURES REACTION'
            total_atoms_products= {k: total_atoms_products.get(k, 0) + cmp_dict.get(k, 0) for k in set(total_atoms_products) | set(cmp_dict)} #adding atoms one by one molecule to the total atoms dictionary

        atoms_difference = {k: total_atoms_reactants.get(k, 0) - total_atoms_products.get(k, 0) for k in set(total_atoms_reactants) | set(total_atoms_products)}

        if 1 in atoms_difference:
            del atoms_difference[1] #removing hydrogen from the account

        for atomic_num, diff in atoms_difference.items():
            if diff != 0:
                return atoms_difference

        return None

    def printAtomsByCompound(self, rxn):
        reactants_with_stoich = [i for i in rxn.split(';') if i.startswith('-')]
        products_with_stoich = [i for i in rxn.split(';') if not i.startswith('-')]

        def expand_compounds_list_with_stoichimetry(comp_list_with_stoich):
            '''
            this funstion expands the list to include molecule as many times as it is mentioned in the stoich coef
            '''
            compounds_list = []
            for cmp_stoich in comp_list_with_stoich:
                stoich = int(cmp_stoich.split(' ')[0])
                comp = cmp_stoich.split(' ')[1]
                compounds_list.extend([comp]*abs(stoich))
            return compounds_list

        reactants = expand_compounds_list_with_stoichimetry(reactants_with_stoich)
        products = expand_compounds_list_with_stoichimetry(products_with_stoich)

        for cmp in reactants:
            print(cmp, self.compound_atoms[cmp])

        for cmp in products:
            print(cmp, self.compound_atoms[cmp])

    def createBalanceFile(self):
        count_unbalanced = 0
        count_unstructured = 0

        # columns of the df_reactions
        #rxnUID,compounds,rxn_stoich_code,balance
        listDF = []
        for index, row in self.df_reactions.iterrows():
            dict_rxn = dict()
            if index%50000 == 0:
                print(index, 'done processing balance')
            rxn = row['rxn_stoich_code']
            dict_rxn['rxnUID'] = row['rxnUID']
            #print(rxn)
            difference_in_balance = self.check_balance(rxn)
            if difference_in_balance == 'EXCLUDE INCLOMPLETE COMPOUND STRUCTURES REACTION':
                count_unstructured += 1
                dict_rxn['balance'] = -1

            elif difference_in_balance:
                count_unbalanced+=1
                dict_rxn['balance'] = 0

            else:
                dict_rxn['balance'] = 1

            listDF.append(dict_rxn)

        f_out = pd.DataFrame(listDF, dtype='int')
        f_out.to_csv(self.data.data_folder+self.data.reaction_network+'/reaction_balance.csv', index=False)

        print('count_unstructured', count_unstructured)
        print('count_unbalanced', count_unbalanced)