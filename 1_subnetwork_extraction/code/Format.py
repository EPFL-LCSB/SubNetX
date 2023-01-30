__author__ = 'anastasia'
import os
import networkx as nx
###############
#  Reactions  #
###############

class Format():

    def __init__(self, data):
        self.data = data
        self.all_compounds = []

    def produceFormattedOutput(self):
        # the main function that produces the input for the optimization
        self.writeReactionFile()
        self.writeCompoundsFile()
        self.createMolfiles()

    ###############
    #  Reactions  #
    ###############

    def writeReactionFile(self):
        with open(self.data.auxilary_output_dir+'/reaction_ids.txt') as f:
            reactions_all = [i.strip() for i in f.readlines()]
        print('number of reactions:', len(reactions_all))

        # writing the reactions in the right format and collecting all compounds
        reaction_file = open(self.data.output_optimization_input_dir+'/reactions.tsv', 'w')
        reaction_file.write('R_PR_UID	R_PR_STOICH\n')

        df_reactions_selected = self.data.df_reactions[self.data.df_reactions['rxnUID'].isin(reactions_all)]
        df_reactions_selected['R_PR_STOICH'] = df_reactions_selected.apply(self.get_equation, axis=1)
        df_reactions_selected.rename(columns={"rxnUID":"R_PR_UID"}, inplace=True)
        df_reactions_selected.drop_duplicates(subset=['R_PR_UID'], inplace=True)
        df_reactions_selected.to_csv(self.data.output_optimization_input_dir+'/reactions.tsv', sep = '\t', columns=['R_PR_UID', 'R_PR_STOICH'], index=False)
        # for rxn in reactions_all:
        #     equation = self.get_equation(rxn)
        #     reaction_file.write('{}\t{}\t{}\n'.format(rxn, equation))
        #     reactions_printed.add(rxn)

    def get_equation(self, row):
        substrates = list()
        products = list()
        for cmp in row['rxn_stoich_code'].split(';'):
            if cmp.startswith('-'):
                substrates.append(cmp.strip('-'))
            else:
                products.append(cmp)
        self.all_compounds.extend([i.split()[1] for i in substrates])
        self.all_compounds.extend([i.split()[1] for i in products])
        equation = ' + '.join(substrates) + ' <=> ' + ' + '.join(products)
        return equation

    ###############
    #  Compounds  #
    ###############

    def writeCompoundsFile(self):

        compounds_all = list(set(self.all_compounds))
        df_compounds_selected = self.data.df_compounds[self.data.df_compounds['cUID'].isin(compounds_all)]
                # columns should be renamed
        # new_df: M_PR_UID	M_PR_NAME	M_PR_INCHIKEY	M_PR_SMILE-CAN-OB	M_PR_FORMULA	M_PR_CHARGE	M_XR_KEGG	M_XR_ALL_OTHER
        # old_df: cUID  COMMON_NAME INCHIKEY SMILES Formula Total_charge

        df_compounds_selected.rename(columns={"cUID": "M_PR_UID",
                                              "COMMON_NAME": "M_PR_NAME",
                                              "INCHIKEY": "M_PR_INCHIKEY",
                                              "SMILES":"M_PR_SMILE-CAN-OB",
                                              "Formula":"M_PR_FORMULA",
                                              "Total_charge": "M_PR_CHARGE"}, inplace = True)
        df_compounds_selected.to_csv(self.data.output_optimization_input_dir+'/compounds.tsv',
                                     columns=["M_PR_UID","M_PR_NAME","M_PR_INCHIKEY","M_PR_SMILE-CAN-OB","M_PR_FORMULA","M_PR_CHARGE"],
                                     sep='\t', index=False)



    def createMolfiles(self):
        if not os.path.exists(self.data.output_optimization_input_dir+'/molfiles'):
            os.mkdir(self.data.output_optimization_input_dir+'/molfiles')
        compounds_file = open(self.data.output_optimization_input_dir+'/compounds.tsv')
        header = compounds_file.readline().split('\t')
        lines = compounds_file.readlines()
        if '' in lines:
            lines.remove('')
        for line in lines:
            smiles = line.split('\t')[header.index('M_PR_SMILE-CAN-OB')]
            os.system('obabel -:"{}" -O {}/molfiles/{}.mol'.format(smiles, self.data.output_optimization_input_dir, line.split('\t')[header.index('M_PR_UID')]))

def createMolfiles(output_optimization_input_dir):
    if not os.path.exists(output_optimization_input_dir+'/molfiles'):
        os.mkdir(output_optimization_input_dir+'/molfiles')
        ready = []
    else:
        ready = [i.replace('.mol','') for i in os.listdir(output_optimization_input_dir)]

    compounds_file = open(output_optimization_input_dir+'/compounds.tsv')
    header = compounds_file.readline().split('\t')
    lines = compounds_file.readlines()
    if '' in lines:
        lines.remove('')
    for line in lines:
        if line.split('\t')[header.index('M_PR_UID')] not in ready:
            smiles = line.split('\t')[header.index('M_PR_SMILE-CAN-OB')]
            os.system('obabel -:"{}" -O {}/molfiles/{}.mol'.format(smiles, output_optimization_input_dir, line.split('\t')[header.index('M_PR_UID')]))

#createMolfiles('/Users/anastasia/Tools/biosubnet/projects/tadalafil/output_optimization_input')