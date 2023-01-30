__author__ = 'anastasia'

from Data import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS
import pandas as pd
import os

class FragmentMatcher():

    def __init__(self):
        self._onPatt = None
        self._offPatts = []

    def AddExclusion(self, sma):
        self._offPatts.append(Chem.MolFromSmarts(sma))

    def InitSML(self, sma):
        self._onPatt = Chem.MolFromSmarts(sma)

    def InitMOL(self, molfile):
        self._onPatt = Chem.MolFromMolFile(molfile)

    def HasMatch(self, mol):
        if self._onPatt is None:
            return 0
        t = mol.HasSubstructMatch(self._onPatt)
        if not t:
            return 0
        else:
            for patt in self._offPatts:
                if mol.HasSubstructMatch(patt):
                    return 0
        return 1

    def GetMatch(self, mol):
        if self._onPatt is None:
            return None
        return mol.GetSubstructMatch(self._onPatt)

    def GetMatches(self, mol, uniquify=1):
        if self._onPatt is None:
            return None
        return mol.GetSubstructMatches(self._onPatt, uniquify=uniquify)

    def GetBond(self, idx):
        if self._onPatt is None:
            return None
        return self._onPatt.GetBondWithIdx(idx)


class compoundFiltering():
    """
    Filter compounds based on predefined substructure
    """
    def __init__(self, data):
        self.data = data

    def getCompoundsWithCarbon(self, list_cmp_ids):
        """
        Filter compound precursors to exclude complex and non-carbon structures
        :param list_cmp_ids:
        :return:
        """
        filtered_list = []
        for cmp in list_cmp_ids:
            if cmp in self.data.dict_smiles.keys():
                smiles = self.data.dict_smiles[cmp]
                carbonsNum = self.numCarbons(smiles)
                if carbonsNum > 0 and carbonsNum < 25:
                    filtered_list.append(cmp)

        return filtered_list

    def numCarbons(self, smiles):
        if type(smiles)==str:
            m = Chem.MolFromSmiles(smiles)
            if m is None: return 0
            numC = 0
            for atom in m.GetAtoms():
                if atom.GetAtomicNum() == 6:
                    numC += 1
            return numC
        return 0

class filterPrecursors(compoundFiltering):

    def __init__(self, data):
        self.data = data
        self.getPrecursorsMol()

    def getPrecursorsMol(self):
        self.precursorMol = dict()
        precursors = self.data.precursorCompounds
        for prec in precursors:
            self.precursorMol[prec]=Chem.MolFromSmiles(self.data.dict_smiles[prec])

    def getMatchingForCompound(self, smiles, fr_matchers):

        sm1 = smiles.replace('[R]','C')
        mol1 = Chem.MolFromSmiles(sm1)
        if mol1:
            opt = []
            for matcher in fr_matchers:
                b1 = matcher.HasMatch(mol1)
                opt.append(b1)
            return any(opt)
        else:
            return False

    def getListMatchingCompounds(self, cmp_name, list_cmp_ids):

        patterns_files = os.listdir('../precursor_patterns/'+cmp_name+'/')
        fr_matchers = []
        for p_file in patterns_files:
            fr_matcher = FragmentMatcher()
            fr_matcher.InitMOL('../precursor_patterns/'+cmp_name+'/'+p_file)
            fr_matchers.append(fr_matcher)

        filtered_list = []
        for cmp in list_cmp_ids:
            if cmp in self.data.dict_smiles.keys():
                smiles = self.data.dict_smiles[cmp]
                if self.getMatchingForCompound(smiles, fr_matchers):
                    filtered_list.append(cmp)
            else: pass

        return filtered_list

    def getMostSimilarCompounds(self, cmp):
        # For methyl-CoA return CoA
        if cmp == '383753545':
            return ['1467865652']

        m = Chem.MolFromSmiles(self.data.dict_smiles[cmp])

        dict_sim_prec = dict()

        for prec, mol_prec in self.precursorMol.items():

            # Return coefficient : num atoms matched / num atoms in precursor
            if mol_prec.GetNumAtoms() > m.GetNumAtoms():
                dict_sim_prec[prec]=0 # we do not need precursor that is more complicated
            else:
                dict_sim_prec[prec]=rdFMCS.FindMCS([mol_prec, m]).numAtoms

        best_match_prec = sorted(dict_sim_prec, key=dict_sim_prec.get, reverse=True)[:self.data.numSimPrecursorsLimit]
        return best_match_prec




