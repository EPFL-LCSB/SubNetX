__author__ = 'anastasia'

import pandas as pd
import os
import shutil
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import hashlib
from itertools import product

class Data():

    def __init__(self):
        """
        read in all the essential data:
        parameters.txt
        network files

        Set up the directories and file paths
        :return: produce data object with all the essential data used through the whole code
        """
        # Handle reading in the project name
        self.current_round_extraction = 0

        if len(sys.argv) == 1:
            print('Please, set project name as the first argument of the script "ex. python3 main.py scoulerine_case1a"')
            print('Alternatively, you can introduce the project name now:')
            self.projectname = input()

        elif len(sys.argv) == 2:
            self.projectname = sys.argv[1]

        else:
            print('Please, set project name as the first argument of the script "ex. python3 main.py scoulerine_case1a"')
            print('Finishing the execution')
            exit()

        # Set folder and file paths
        self.data_folder = '../data/'

        self.excludelist_compounds = '../projects/'+self.projectname+'/excludelists/compounds.txt'
        self.toxic_compounds_file = '../projects/'+self.projectname+'/excludelists/toxic_compounds.txt'
        self.mammal_cofactor_file = '../projects/'+self.projectname+'/excludelists/mammal_cofactors.txt'

        self.excludelist_reactions = '../projects/'+self.projectname+'/excludelists/reactions.txt'

        self.compound_limits_file = '../projects/'+self.projectname+'/compound_parameters.csv'
        self.NICEpathways = '../projects/'+self.projectname+'/auxilary_output/NICEpathways_linear.txt'

        # Setting the defaults in case the individual parameter files were not provided
        if os.path.exists( '../projects/'+self.projectname+'/parameters.txt'):
            self.param_file = '../projects/'+self.projectname+'/parameters.txt'
        else:
            self.param_file = '../defaults/parameters.txt'

        if not os.path.exists(self.compound_limits_file):
            self.compound_limits_file = '../defaults/compound_parameters.csv'

        if not os.path.exists(self.excludelist_compounds):
            self.excludelist_compounds = "../defaults/excludelists/compounds.txt"

        if not os.path.exists(self.toxic_compounds_file):
            self.toxic_compounds_file = "../defaults/excludelists/toxic_compounds.txt"

        if not os.path.exists(self.mammal_cofactor_file):
            self.mammal_cofactor_file = "../defaults/excludelists/mammal_cofactors.txt"

        if not os.path.exists(self.excludelist_reactions):
            self.excludelist_reactions = "../defaults/excludelists/reactions.txt"


        with open(self.excludelist_compounds) as f:
           self.excludecompounds = [i.strip() for i in f.readlines() if i!= '\n']

        with open(self.toxic_compounds_file) as f:
           self.toxiccompounds = [i.strip().split()[0] for i in f.readlines() if i!= '\n']

        with open(self.mammal_cofactor_file) as f:
            self.mammalcofactors = [i.strip().split()[0] for i in f.readlines() if i!= '\n']

        with open(self.excludelist_reactions) as f:
           self.excludereactions = [i.strip() for i in f.readlines() if i!= '\n']


        # create lists to collect boundaries and branching points
        self.boundary_no_path = [] #
        self.no_branching_points = []
        self.branching_points_1 = []
        self.branching_points_2plus = []
        self.precursorCompounds = []
        self.boundaries_total= []

        # set the default parameters
        self.readParametersFile("../defaults/parameters.txt")
        # reset for the specific project parameters
        self.readParametersFile("../projects/"+self.projectname+"/parameters.txt")

        # set the metabolites of the model
        self.getModelMetabolites()

        # Set paths to .csv files of the reaction network
        self.setNetworksPath()

        # Read in the networks (e.g. ARBRE)
        self.readInitialDataframes()

        # if auxilary networks required append them to the main dataframes
        if self.use_auxilary_network:
            self.readAuxilaryDataframes()
            self.processAuxilaryDataframes()

        self.processDataframes()

        # create directories for the output
        self.createDirectories()

    ####################################################################################################################
    #                                    set parameters and paths
    ####################################################################################################################

    def readParametersFile(self, filename):
        """
        Set the parameters of the execution
        :return:
        """
        dict_param = dict()
        with open(filename) as f:
            for line in f:
                if not line.startswith('#') and line.strip() != '':
                    dict_param[line.split('|')[0]] = line.split('|')[1].strip()
        self.model_organism = dict_param['model_organism']
        if dict_param['main_precursor'] == 'all': self.main_precursor = dict_param['main_precursor']
        else: self.main_precursor = dict_param['main_precursor']
        self.main_target = dict_param['main_target']
        self.num_shortest_pathways = int(dict_param['num_shortest_pathways'])
        if 'minplus' in dict_param: self.minplus = int(dict_param['minplus'])
        self.num_pathways_to_model = int(dict_param['num_pathways_to_model'])
        self.run_expansion = [True if dict_param['run_expansion']=='1' else False][0]
        self.filter_precursor_structure = [True if dict_param['filter_precursor_structure']=='1' else False][0]
        self.lowest_atom_conservation_threshold=float(dict_param['lowest_atom_conservation_threshold'])
        if 'structure_based_pairs' in dict_param: self.structure_based_pairs = int(dict_param['structure_based_pairs'])

        if 'prefer_known' in dict_param:self.prefer_known = int(dict_param['prefer_known'])
        if 'use_exponential_transformation' in dict_param: self.use_exponential_transformation = int(dict_param['use_exponential_transformation'])
        # possible distance types: dist, dist_exp, dist_known, dist_exp_known
        if self.prefer_known == 1 and self.use_exponential_transformation == 1:
            self.distance_transformation = 'dist_exp_known' # known exponential distance
        elif self.prefer_known == 1 and self.use_exponential_transformation == 0:
            self.distance_transformation = 'dist_known' # known distance
        elif self.prefer_known == 0 and self.use_exponential_transformation == 1:
            self.distance_transformation = 'dist_exp' # exponential distance
        elif self.prefer_known == 0 and self.use_exponential_transformation == 0:
            self.distance_transformation = 'dist' # standard distance

        self.boundaries_alternatives_num=int(dict_param['boundaries_alternatives_num'])
        self.numSimPrecursorsLimit=int(dict_param['numSimPrecursorsLimit'])

        self.reaction_network = dict_param['reaction_network']
        if 'use_auxilary_network' in dict_param: self.use_auxilary_network = [True if dict_param['use_auxilary_network']=='1' else False][0]

        self.graph_file = self.data_folder+self.reaction_network+'/'+self.reaction_network+'.gpickle'

    def createDirectories(self):
        # create directories for figures, stats and output networks
        if not os.path.exists('../projects/'+self.projectname):
            os.mkdir('../projects/'+self.projectname)

        # Collecting the statistics of the network extraction: generations, number of compounds, number of reactions etc
        self.stats_dir = '../projects/'+self.projectname+'/stats/'
        if not os.path.exists(self.stats_dir):
            os.mkdir(self.stats_dir)

        # the output of the network extraction that will be passed to the optimization is in this folder
        self.output_optimization_input_dir = '../projects/'+self.projectname+'/output_optimization_input/'
        if not os.path.exists(self.output_optimization_input_dir):
            os.mkdir(self.output_optimization_input_dir)

        # visualizations are stored in this folder
        self.fig_dir = '../projects/'+self.projectname+'/figures/'
        if not os.path.exists(self.fig_dir):
            os.mkdir(self.fig_dir)

        # intermediate results that can be used for detailed analysis are stored here
        self.auxilary_output_dir = '../projects/'+self.projectname+'/auxilary_output/'
        if not os.path.exists(self.auxilary_output_dir):
            os.mkdir(self.auxilary_output_dir)
        self.reaction_ids_file = self.auxilary_output_dir + "reaction_ids.txt"

        # the compounds that are identified as boundaries are stored here
        self.boundaries_dir = '../projects/'+self.projectname+'/auxilary_output/boundaries'
        if not os.path.exists(self.boundaries_dir):
            os.mkdir(self.boundaries_dir)
        else:
            shutil.rmtree(self.boundaries_dir)
            os.mkdir(self.boundaries_dir)

    def setNetworksPath(self):
        self.network_file = self.data_folder+self.reaction_network+'/network.csv'
        self.reactions_pairs_file = self.data_folder+self.reaction_network+'/reactions_pairs.csv'
        self.reactions_file = self.data_folder+self.reaction_network+'/reactions.csv'
        self.comp_file = self.data_folder+self.reaction_network+'/compounds.csv'
        self.reaction_balance = self.data_folder+self.reaction_network+'/reaction_balance.csv'
        # auxilary networks - not always present
        self.auxilary_reactions_file = '../projects/'+self.projectname+'/auxilary_network/reactions.csv'
        self.auxilary_network_file = '../projects/'+self.projectname+'/auxilary_network/network.csv'

    ####################################################################################################################
    #                       get model and initial pathway metabolites
    ####################################################################################################################

    def getModelMetabolites(self):
        metabolitesFile = self.data_folder+'organisms_metabolites_annotated/'+self.model_organism+'.tsv'
        df_model_metabolites = pd.read_csv(metabolitesFile, sep = '\t')
        df_model_metabolites.dropna(subset=['metaboliteLCSBID'], inplace=True)
        self.modelMetabolites = set(df_model_metabolites['metaboliteLCSBID'].to_list()) - {'None'}
        print('Model metabolites total: ', len(self.modelMetabolites))

    def getInitialPathways(self):
        init_pathways = []
        with open(self.NICEpathways) as f:
            for line in f:
                init_pathways.append(line.split('->'))
        return init_pathways

    ####################################################################################################################
    #                                        Reading in the dataframes
    ####################################################################################################################

    def readInitialDataframes(self):
        """
        Reading the data into the pandas dataframes
        :return:
        """
        self.df_network = pd.read_csv(self.network_file)
        self.df_network = self.df_network.astype({'Instance ID of pair':str,'UID of pair':str,'source':str, 'target':str})

        self.df_reactions = pd.read_csv(self.reactions_pairs_file)
        self.df_reactions = self.df_reactions.astype({'Instance ID of pair':str,'UID of pair':str,'rxnUID':str})

        self.df_reactions_comp = pd.read_csv(self.reactions_file)
        self.df_reactions_comp = self.df_reactions_comp.astype({'rxnUID':str})

        self.df_compounds = pd.read_csv(self.comp_file, low_memory=False)
        self.df_compounds = self.df_compounds.astype({'cUID': str})

        self.df_reaction_balance = pd.read_csv(self.reaction_balance)
        self.df_reaction_balance = self.df_reaction_balance.astype({'rxnUID':str})

        self.df_compound_properties_limits = pd.read_csv(self.compound_limits_file)

    def readAuxilaryDataframes(self):
        """
        Reading the data into the pandas dataframes
        :return:
        """
        self.df_auxilary_network = pd.read_csv(self.auxilary_network_file)
        self.df_auxilary_reactions = pd.read_csv(self.auxilary_reactions_file)

    ####################################################################################################################
    #                                           Processing dataframes
    ####################################################################################################################

    def removeReactionsCompounds(self, mode='toxic'):
        """
        Remove tox compounds for all, mammal compounds for yeast and ecoli
        :param mode:
        :return:
        """
        if mode == 'toxic':
            self.df_reactions_comp['to_remove'] = self.df_reactions_comp.apply(self.markToxicReactions, axis = 1)
        elif mode == 'mammal':
            self.df_reactions_comp['to_remove'] = self.df_reactions_comp.apply(self.markMammalCofactorReactions, axis = 1)
        self.df_reactions_comp = self.df_reactions_comp[self.df_reactions_comp['to_remove']==0]

    def markToxicReactions(self, row):
        if any(i in row['compounds'].split(';') for i in self.toxiccompounds):
            return 1
        return 0

    def markMammalCofactorReactions(self, row):
        if any(i in row['compounds'].split(';') for i in self.mammalcofactors):
            return 1
        return 0

    def processDataframes(self):
        """
        Exclude unbalanced reactions
        Add pairs to reactions
        Add scores to reactions
        :return:
        """
        # self.removeReactionsCompounds(mode='toxic')
        if self.model_organism == 'ecoli' or self.model_organism == 'yeast':
            self.removeReactionsCompounds(mode='mammal')

         # filter out the reactions that are in excludelists/reactions.txt
        self.df_reactions = self.df_reactions[~self.df_reactions.rxnUID.isin(self.excludereactions)]
        
        # keep only balanced reactions
        self.df_reaction_balance = self.df_reaction_balance[self.df_reaction_balance['balance']==1]
        self.df_reactions = self.df_reactions.merge(self.df_network, on="UID of pair", how='inner')
        self.df_reactions = self.df_reactions.merge(self.df_reaction_balance, on='rxnUID', how='inner')
        self.df_reactions = self.df_reactions.merge(self.df_reactions_comp, on='rxnUID', how='inner')
        
        edges_uids_keep = set(self.df_reactions['UID of pair'].to_list())
        self.df_network = self.df_network[self.df_network['UID of pair'].isin(edges_uids_keep)]

        # exclude / keep pairs with structure based similarity score
        self.df_network = self.df_network[self.df_network['structure_based']<=self.structure_based_pairs]

        # get dictionary with smiles for each compound
        self.dict_smiles = dict(zip(self.df_compounds.cUID, self.df_compounds.SMILES))

    def processAuxilaryDataframes(self):
        """
        Transform the generic input format from reactions and network into the same format as the total input and
        append it to the total network
        :return: append new compounds to df compounds
        """
        #
        self.newCompounds()
        self.newReactions()
        self.newReactionsBalance()
        self.newNetwork()
        self.newReactionsPairs()

    def newReactions(self):
        # rxnUID,compounds,rxn_stoich_code,if_known
        if not os.path.exists('../projects/'+self.projectname+'/auxilary_network/reactions_aux.csv'):
            list_df_out=[]

            for rxn_inchikey in self.df_auxilary_reactions['REACTION_INCHIKEY'].to_list():

                inchikey_compounds = []
                stoich_compounds = []

                inchikey_compounds.extend([i.split()[1] for i in rxn_inchikey.split(' <=> ')[0].split(' + ')]) # add reactants
                stoich_compounds.extend([int(i.split()[0])*-1 for i in rxn_inchikey.split(' <=> ')[0].split(' + ')]) # add reactants stoich

                inchikey_compounds.extend([i.split()[1] for i in rxn_inchikey.split(' <=> ')[1].split(' + ')]) # add products
                stoich_compounds.extend([i.split()[0] for i in rxn_inchikey.split(' <=> ')[1].split(' + ')]) # add products stoich

                lcsb_ids = []
                for inchikey in inchikey_compounds:
                    lcsb_ids.append(self.getLCSBidFromInchikey(inchikey))

                list_df_out.append({
                    "rxnUID": hashlib.sha256(rxn_inchikey.encode()).hexdigest(),
                    "compounds": ';'.join(str(i) for i in set(lcsb_ids)),
                    "rxn_stoich_code": ';'.join([str(stoich_compounds[i])+' '+str(lcsb_ids[i]) for i in range(len(lcsb_ids))]),
                    "if_known": 0
                })

            df_react_aux = pd.DataFrame(list_df_out)

            df_react_aux.to_csv('../projects/'+self.projectname+'/auxilary_network/reactions_aux.csv', index=False)

        else:
            df_react_aux = pd.read_csv('../projects/'+self.projectname+'/auxilary_network/reactions_aux.csv')
        self.df_reactions_comp = pd.concat([self.df_reactions_comp, df_react_aux])

    def newReactionsPairs(self):
        # Instance ID of pair,UID of pair,rxnUID,if_known_check,max_bridgit
        if not os.path.exists('../projects/'+self.projectname+'/auxilary_network/reactions_pairs_aux.csv'):
            list_df_out=[]

            count = 1

            for rxn_inchikey in self.df_auxilary_reactions['REACTION_INCHIKEY'].to_list():

                rxnUID = hashlib.sha256(rxn_inchikey.encode()).hexdigest()

                inchikey_reactants = [i.split()[1] for i in rxn_inchikey.split(' <=> ')[0].split(' + ')]
                inchikey_products = [i.split()[1] for i in rxn_inchikey.split(' <=> ')[1].split(' + ')]

                all_pairs = product(inchikey_reactants, inchikey_products)

                for pair in all_pairs:
                    pair_hash = self.getPairInchikeyHash(list(pair))
                    if pair_hash in self.df_auxilary_network['UID of pair'].to_list():
                        list_df_out.append({
                            'Instance ID of pair':count,
                            'UID of pair':pair_hash,
                            'rxnUID':rxnUID,
                            'if_known_check':0,
                            'max_bridgit':1
                            })
                    count += 1
            df_react_aux = pd.DataFrame(list_df_out)
            df_react_aux.to_csv('../projects/'+self.projectname+'/auxilary_network/reactions_pairs_aux.csv', index=False)
        else:
            df_react_aux = pd.read_csv('../projects/'+self.projectname+'/auxilary_network/reactions_pairs_aux.csv')
        self.df_reactions = pd.concat([self.df_reactions, df_react_aux])

    def getLCSBidFromInchikeyRow(self, row, col_name):
        return self.getLCSBidFromInchikey(row[col_name])

    def getLCSBidFromInchikey(self, inchikey):
        df_comp_inchikey = self.df_compounds[self.df_compounds['INCHIKEY']==inchikey]
        if len(df_comp_inchikey)>0:
            return df_comp_inchikey['cUID'].iloc[0]
        return inchikey

    def newReactionsBalance(self):
        # rxnUID,balance
        if not os.path.exists('../projects/'+self.projectname+'/auxilary_network/reactions_balance_aux.csv'):
            list_df_out=[]

            for rxn_inchikey in self.df_auxilary_reactions['REACTION_INCHIKEY'].to_list():
                list_df_out.append({
                    "rxnUID": hashlib.sha256(rxn_inchikey.encode()).hexdigest(),
                    "balance": 1
                })
            df_react_aux = pd.DataFrame(list_df_out)
            df_react_aux.to_csv('../projects/'+self.projectname+'/auxilary_network/reactions_balance_aux.csv', index=False)
        else:
            df_react_aux = pd.read_csv('../projects/'+self.projectname+'/auxilary_network/reactions_balance_aux.csv')

        self.df_reaction_balance = pd.concat([self.df_reaction_balance, df_react_aux])

    def newNetwork(self):
        # Instance ID of pair,UID of pair,score,source,target,structure_based,dist,dist_known,dist_exp,dist_exp_known,known_reaction,max_bridgit
        #if not os.path.exists('../projects/'+self.projectname+'/auxilary_network/network_aux.csv'):
        self.df_auxilary_network['Instance ID of pair'] = ['aux'+str(i+1) for i in range(len(self.df_auxilary_network))]
        self.df_auxilary_network['UID of pair'] = self.df_auxilary_network.apply(self.getPairInchikeyHashRow, axis=1)
        self.df_auxilary_network['source'] = self.df_auxilary_network.apply(self.getLCSBidFromInchikeyRow, args=('CMP1',), axis=1)
        self.df_auxilary_network['target'] = self.df_auxilary_network.apply(self.getLCSBidFromInchikeyRow, args=('CMP2',), axis=1)
        self.df_auxilary_network['structure_based'] = 0
        self.df_auxilary_network['dist'] = self.df_auxilary_network.apply(self.getDistAux, axis=1)
        self.df_auxilary_network['dist_known'] = self.df_auxilary_network.apply(self.getDistAux, axis=1)
        self.df_auxilary_network['dist_exp'] = self.df_auxilary_network.apply(self.getDistAux, axis=1)
        self.df_auxilary_network['dist_exp_known'] = self.df_auxilary_network.apply(self.getDistAux, axis=1)
        self.df_auxilary_network['known_reaction']=0
        self.df_auxilary_network['max_bridgit']=1
        self.df_auxilary_network.rename(columns={'SCORE':'score'}, inplace=True)
        self.df_auxilary_network.to_csv('../projects/'+self.projectname+'/auxilary_network/network_aux.csv',
                                        columns = [
                                            'Instance ID of pair',
                                            'UID of pair',
                                            'score',
                                            'source',
                                            'target',
                                            'structure_based',
                                            'dist',
                                            'dist_known',
                                            'dist_exp',
                                            'dist_exp_known',
                                            'known_reaction',
                                            'max_bridgit'
                                        ], index=False)
        df_network_aux = pd.read_csv(('../projects/'+self.projectname+'/auxilary_network/network_aux.csv'))
        self.df_network = pd.concat([self.df_network, df_network_aux])

    def getDistAux(self, row):
        return 1/row['SCORE']

    def getPairInchikeyHashRow(self, row):
        cmps = [row['CMP1'], row['CMP2']]
        return self.getPairInchikeyHash(cmps)

    def getPairInchikeyHash(self, cmps):
        cmps.sort()
        return hashlib.sha256(''.join(cmps).encode()).hexdigest()

    def newCompounds(self):# collect new compounds
        if not os.path.exists('../projects/'+self.projectname+'/auxilary_network/compounds_aux.csv'):
            smiles_compounds = []
            for rxn_smiles in self.df_auxilary_reactions['REACTION_SMILES'].to_list():
                smiles_compounds.extend(rxn_smiles.split('>>')[0].split('.')) # add reactants
                smiles_compounds.extend(rxn_smiles.split('>>')[1].split('.')) # add products
            smiles_compounds = set(smiles_compounds)
            list_df_out=[]
            dict_smiles_inchikey = dict()
            for smiles in smiles_compounds:
                dict_smiles_inchikey[smiles] = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
                if not dict_smiles_inchikey[smiles] in self.df_compounds['INCHIKEY']:
                    list_df_out.append(self.processNewCompound(dict_smiles_inchikey, smiles))
            df_comp_aux = pd.DataFrame(list_df_out)
            df_comp_aux.to_csv('../projects/'+self.projectname+'/auxilary_network/compounds_aux.csv', index=False)
        else:
            df_comp_aux = pd.read_csv('../projects/'+self.projectname+'/auxilary_network/compounds_aux.csv')
        self.df_compounds = pd.concat([self.df_compounds, df_comp_aux])

    def processNewCompound(self, dict_smiles_inchikey, smiles):
        # generate dictionary with new compound description
        # output file header:
        # cUID,NUM_PATENTIDS,NUM_PUBMEDIDS,NUM_TOTAL,COMMON_NAME,OTHER_NAMES,OLD_ENTRIES,Formula,
        # if_positive_charge,if_negative_charge,if_triple_bond,smiles_length,Num_atoms,MW,Num_Rings,Num_rotatable_bonds,
        # Total_charge,benzene,pyridine,pyrimidine,pyridazine,imidazole,imidazole_cation,pyrrole,pyrazine,furan,pyrazole,
        # thiophene,isoxazole,isoxazole_cation,oxazole,thiazole,thiazole_cation,tetrazole_1234,tetrazole_1235,triazole_134,
        # triazole_123,triazole_124,thiadiazole_134,thiadiazole_cation_123,undefined_ring_type,total_negative_charge,
        # total_positive_charge,H,Cu,Co,Mg,I,F,C,Cl,B,Se,As,Na,Ni,Mn,Ca,Pb,Br,N,Cr,Cd,Mo,K,Ag,Hg,Fe,O,P,S,Li,Te,Si,SMILES,
        # NucleotideRTransferMolecule,Polymer,undefined_structure,INCHIKEY,if_CoA_compound,LD50_mg_kg
        dict_comp = dict()
        for col in self.df_compounds.columns:
            dict_comp[col]=0 # make dummy filling to each column to later
        dict_comp['cUID']=dict_smiles_inchikey[smiles].split('-')[0]
        dict_comp['Formula']=rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))
        dict_comp['MW']=rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(smiles))
        dict_comp['Num_atoms']=Chem.MolFromSmiles(smiles).GetNumAtoms()
        dict_comp['INCHIKEY']=dict_smiles_inchikey[smiles].split('-')[0]
        dict_comp['SMILES']=smiles
        dict_comp['smiles_length']=len(smiles)
        dict_comp['LD50_mg_kg']=10000 # dummy toxicity estimate that fits into the standard range
        return dict_comp

    ####################################################################################################################
    #                                          Writing the output
    ####################################################################################################################

    def writeReactionsAll(self, G):
        reactions_all = []
        for edge in G.edges():
            reactions_all.extend(list(G[edge[0]][edge[1]]['reactions']))

        reactions_all = set(reactions_all)
        with open(self.reaction_ids_file, 'w') as reactions_file:
            for reaction in reactions_all:
                reactions_file.write('{}\n'.format(reaction))

    def writeOutputBoundaryInfoList(self):
        with open(self.boundaries_dir+'/boundary_no_path.txt', 'w') as w:
            w.write('\n'.join([str(i) for i in self.boundary_no_path]))

        with open(self.boundaries_dir+'/no_branching_points.txt', 'w') as w:
            for br_p in self.no_branching_points:
                w.write('{}\n'.format(br_p))

        with open(self.boundaries_dir+'/branching_points_1_comp.txt', 'w') as w:
            for br_p in self.branching_points_1:
                w.write('{}\n'.format(br_p))

        with open(self.boundaries_dir+'/branching_points_2plus_comp.txt', 'w') as w:
            for br_p in self.branching_points_2plus:
                w.write('{}\n'.format(br_p))