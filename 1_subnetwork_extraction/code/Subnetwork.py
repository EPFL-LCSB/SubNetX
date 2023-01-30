__author__ = 'anastasia'

import networkx as nx
from itertools import islice
import numpy as np
import pandas as pd
import os
from datetime import datetime
import matplotlib.pyplot as plt
from Compound import *
pd.options.mode.chained_assignment = None

class Graph():

    def __init__(self, data):
        self.G = nx.Graph()
        self.data = data

    ########  Functions for transformations of the graph preliminary to pathway search #######

    def formG(self):
        """
        Form networkx graph object from pandas df
        """
        graph_file = self.data.graph_file
        # Check if the graph file exists:
        # ATTENTION : if compound filtering or any other graph parameters changed, delete the graph file for it will not
        # be updated
        if os.path.exists(graph_file):
            self.G =  self.getGraph(graph_file)
        elif os.path.exists(self.data.network_file
                 and os.path.exists(self.data.reactions_file)
                 and os.path.exists(self.data.reaction_balance)):
            self.createGraphFile(graph_file)
        else:
            print('PARAMETER ERROR: Defined in the parameters graph file was not found')
            exit()

        print('Graph before filtering compounds: ', len(self.G))

        compoundsPassFilter = self.getFittingCompounds()

        self.G = self.getSubgraphCompounds(compoundsPassFilter)

        print('Graph after filtering compounds: ', len(self.G))

        print("Graph loaded")

    def createGraphFile(self, graph_file):
        """
        :param graph_file:
        :return:
        """
        self.G = nx.Graph()

        print("Load network...")
        print(len(self.data.df_reactions))
        df_network = self.data.df_reactions[self.data.df_reactions['score'] >= self.data.lowest_atom_conservation_threshold]

        for index, row in df_network.iterrows():
            self.G.add_edges_from([(str(row['source']), str(row['target']),
                                        {'id': row['UID of pair'], 'car': row['score'],
                                         'dist': row['dist'], 'dist_known':row['dist_known'],
                                         'dist_exp': row['dist_exp'], 'dist_exp_known':row['dist_exp_known']})])
        print("Write network for future iterations...")
        nx.write_gpickle(self.G, graph_file)

    def getFittingCompounds(self):
        """
        This function truncates the list of compounds to the compounds with the properties defined in the compound_parameters.csv file
        :return: list of LCSB IDs of the compounds that fit the parameter properties
        """
        df_compound_limits = self.data.df_compound_properties_limits
        query_list = []
        for index, row in df_compound_limits.iterrows():
            # e.g.'(Num_atoms>1&Num_atoms<1000) and (MW>0&MW<1000)'
            query_list.append(row['Code']+'>='+str(row['Min'])+'&'+row['Code']+'<='+str(row['Max']))

        query = '&'.join(query_list)

        df_filtered_compounds = self.data.df_compounds.query(query, engine="python")

        cUIDs = df_filtered_compounds['cUID'].to_list()

        # excude undesired and toxic compounds
        cUIDs = list(set(cUIDs)-set(self.data.excludecompounds)-set(self.data.toxiccompounds))

        return cUIDs

    def getSubgraphCompounds(self, cUIDs):
        H = self.G.subgraph(cUIDs)
        return H


    def getGraph(self, filename):
        """
        read graph from gpickle file
        :param filename: name of the .gpickle file that was created in the precious iterations of running the algorithm
        :return: networkx graph object
        """
        G = nx.read_gpickle(filename)
        return G

    #############            Functions for pathway search within Graph             ###############

    def setNumberOfPathways(self, num_pathways):
        """
        This function sets number of shortest pathways that will be found within the initial graph
        :param num_pathways:
        """
        self.number_of_pathways = num_pathways

    def setDistanceType(self, distance_type):
        self.distance_type = distance_type

    def findPrecursorCompoundsInGraph(self):
        """
        Identify as precursors only those compounds that are in the graph
        :return:
        """
        self.precursorCompounds = self.data.precursorCompounds.intersection(set(self.G.nodes()))
        print("Number of model metabolites in the network:", len(self.precursorCompounds))

    def findInitialPathways(self, precursor, target, num_pathways):
        """
        Use linear pathway search to find the core set of pathways
        :param precursor: precursor compound id
        :param target: target compound id
        :param num_pathways: total number of linear pathways to find
        :return:
        """
        allPathways = self.k_shortest_paths(precursor, target, num_pathways, self.distance_type)
        return allPathways

    def findShortestPathwaysToModel(self, target):
        #collecting all shortest pathways to a target from each precursor in the model
        allPathways = []
        fiP = filterPrecursors(self.data) # initiate an object to filter precursors to find the most similar one
        precs = fiP.getMostSimilarCompounds(target)
        for precursor in  precs:
            allPathways.extend(self.k_shortest_paths(precursor, target, self.number_of_pathways, self.distance_type))

        if not allPathways:
            # add this compound to the list of compounds for which no pw to any metabolite exists in the network
            self.data.boundary_no_path.append(target)
            return allPathways

        allPathways = self.getMinimalPathways(allPathways)

        return allPathways

    def getMinimalPathways(self, pathwaysOneTarget):
        if not pathwaysOneTarget:
            return []
        min_length = min(len(path) for path in pathwaysOneTarget)
        length_limit = self.data.minplus + min_length
        totalPathways = [path for path in pathwaysOneTarget if len(path)<=length_limit]
        return totalPathways

    def k_shortest_paths(self, source, target, k, weight):
        if target == '383753545':
            # for methyl CoA only return 1 (1-step hydrolysis) pathway towards CoA
            return list(islice(nx.shortest_simple_paths(self.G, source, target, weight=weight), 1))
        try:
            return list(islice(nx.shortest_simple_paths(self.G, source, target, weight=weight), k))
        except Exception as e:
            print('Exception happened while searching for k_shortest_path (Subnetwork.py)')
            print(e)
            return []

    def writeGraphGpickle(self, file):
        nx.write_gpickle(self.G, file)

class Subnetwork(Graph):

    def addPathwayFromPathwaySearchOutputToSubnetwork(self, pathway_list):
        intermediates = [str(i) for i in pathway_list]
        pairs = self.getPairsIDs(intermediates)
        for pair in pairs:
            cmp1, cmp2 = pair[0], pair[1]
            try:
                self.G.add_edge(cmp1, cmp2)
            except Exception as e:
                print(e)

    def getPairsToAnnotate(self):
        self.pairs_to_annotate = []
        for edge in self.G.edges:
            if 'reactions' not in self.G[edge[0]][edge[1]].keys():
                self.pairs_to_annotate.append((edge[0], edge[1]))

    def annotateNewPairsWithReactionsAndBoundary(self):
        """
        This function is imortant as it calculates what are the bounaries that we care about
        :return:
        """
        for pair in self.pairs_to_annotate:
            all_boundary=[] # create a list to collect all boundary that will be collected for this pair according to parameters
            rxnUIDs=[] # create a list of reactions to collect all reactions that are associated to the boundaries
            flag_branching = True # flag identifying what is the minimal branching of this point
            compound1 = pair[0]
            compound2 = pair[1]

            # select only the data related to the pair that we are annotating
            df_pair = self.data.df_reactions.query(
                'source==@compound1&target==@compound2|source==@compound2&target==@compound1')

            # sanity check - if all pairs are assigned at least 1  reaction
            rxn_alternatives = len(df_pair['rxnUID'].to_list()) # total number of reaction alternatives
            if rxn_alternatives==0:
                print('no reactions for pair', pair)
                return False

            # calculating boundary compounds set
            # here are all boundaries - later is boundary filtering based on branching points
            df_pair['boundary'] = df_pair.apply(self.getBoundaryPerReaction, axis=1, args=(compound1, compound2,))
            df_pair['numBoundary'] =  df_pair.apply(self.getNumBoundaryPerReaction, axis=1)

            alternatives_count_left = self.data.boundaries_alternatives_num

            # if we have an option to take a reaction without boundaries, we take it
            df_pair_no_boundary = df_pair[df_pair['numBoundary']==0] # zero boundaries = no branching
            if len(df_pair_no_boundary) != 0:
                df_slice = df_pair_no_boundary.iloc[0:alternatives_count_left+1, :] # take slice according to the defined number of rxn alternatives
                rxnUIDs.extend(df_slice['rxnUID'].to_list()) # add reactions that are in this slice
                all_boundary.extend([i for i in df_slice['boundary'].to_list() if i!='']) # add boundaries for this reaction to all boundaries list
                alternatives_count_left -= len(df_slice) #
                self.data.no_branching_points.append(pair)
                flag_branching = False


            # otherwise, we append reactions with only 1 boundary
            if flag_branching == True:
                df_pair_1_boundary = df_pair[df_pair['numBoundary']==1]
                if len(df_pair_1_boundary) != 0 and alternatives_count_left > 0:
                    df_slice = df_pair_1_boundary.iloc[0:alternatives_count_left+1, :]
                    rxnUIDs.extend(df_slice['rxnUID'].to_list())
                    all_boundary.extend([i for i in df_slice['boundary'].to_list() if i!=''])
                    alternatives_count_left -= len(df_slice)
                    self.data.branching_points_1.append(pair)
                    flag_branching = False


            # otherwise, we have to take reactions with more that 1 boundary
            if flag_branching == True:
                df_pair_2plus_boundary = df_pair[df_pair['numBoundary']>1]
                if len(df_pair_2plus_boundary) != 0  and alternatives_count_left > 0:
                    df_slice = df_pair_2plus_boundary.iloc[0:alternatives_count_left+1, :]
                    rxnUIDs.extend(df_slice['rxnUID'].to_list())
                    for compound_set in df_slice['boundary'].to_list():
                        if set(compound_set.split(',')) != set(['1467880389','1469319694']): # cocaine to homocaine cannot be cofactors! pairUID 2603588252
                            all_boundary.extend([i for i in compound_set.split(',')])
                    alternatives_count_left -= len(df_slice)
                    self.data.branching_points_2plus.append(pair)

            # saving all reactions and boundary to the graph
            self.G[pair[0]][pair[1]]['reactions'] = set(rxnUIDs)
            self.G[pair[0]][pair[1]]['all_rxns_boundary'] = set(all_boundary)

        return True

    def getBoundaryPerReaction(self, row, compound1, compound2):
        """
        :param row: row of the dataframe that contains boundaries per reactant-product pair per reaction
        :param compound1: compound of the reactant-product pair
        :param compound2: compound of the reactant-product pair
        :return:
        """
        return ','.join(
            [str(i) for i in list(
                set(row['compounds'].split(';')) - {compound1, compound2} - self.data.modelMetabolites
            )]
        )

    def getNumBoundaryPerReaction(self, row):
        return(len(row['boundary'].split(',')))

    def getPairsIDs(self, intermediates):
        pairs = []
        for i in range(len(intermediates)-1):
            pairs.append((intermediates[i], intermediates[i+1]))
        return pairs

    def identifyOuterBoundaries(self):
        all_boundaries = []
        for edge in self.G.edges:
            all_boundaries.extend(self.G[edge[0]][edge[1]]['all_rxns_boundary'])

        # 1468093962 : H2 remove as it is gas
        # 1468141273 : a [4Fe-4S] iron-sulfur cluster
        # 1468141265 : a [4Fe-4S] iron-sulfur cluster
        # 1467879821 : sodium cation
        # 1468034880 : cobalt
        # 1468083953 : Cu2+
        # 1467881253 : Hydrogen selenide
        # 1467873656 : N2 nitrogen
        # 1468339327 : Co
        # 1467869449 : Mg2+
        return set(all_boundaries) - set(self.G.nodes) - set(self.data.modelMetabolites) \
               - {'1468093962', '1468141273', '1468141265', '1467879821', '1468034880', '1468083953', '1467881253',
                  '1467873656', '1468339327', '1467869449'}

    # ===============================================================================
    #
    # Calculate statistics for the resulting subnetwork ---------->

    def calculateStatistics(self):
        dateAndTime = datetime.now().isoformat(timespec='minutes')
        graph_stats_file_name = 'converged_graph_stats.csv'
        outfile = open(self.data.stats_dir+'/'+graph_stats_file_name, 'w')

        outfile.write('num_nodes,%d\n'%(self.G.number_of_nodes()))
        outfile.write('num_edges,%d\n' % (self.G.number_of_edges()))
        outfile.write('num_nodes_connected,%d\n' % ((self.G.number_of_nodes() - len(list(nx.isolates(self.G))))))
        outfile.write('num_isolated_islands,%d\n' % (len(list(self.G.subgraph(c) for c in nx.connected_components(self.G)))))

        # Statistics on main / giant component
        if len(self.G) > 0:
            main_component_G = max([self.G.subgraph(c) for c in nx.connected_components(self.G)], key=len)
            outfile.write('num_nodes_biggest_island,%d\n' % (main_component_G.number_of_nodes()))
            outfile.write('num_edges_biggest_island,%d\n' % (main_component_G.number_of_edges()))
            outfile.write('percent_nodes_biggest_island,%s\n' % (str(round((main_component_G.number_of_nodes() / self.G.number_of_nodes() * 100),2))))
            outfile.write('percent_edges_biggest_island,%s\n' % (str(round((main_component_G.number_of_edges() / self.G.number_of_edges() * 100),2))))
            ifTargetInMain = self.data.main_target in main_component_G.nodes
            outfile.write('if_target_in_biggest_island,{}\n'.format(ifTargetInMain))
            if self.data.main_precursor != 'all':
                ifPrecursorInMain = self.data.main_precursor in main_component_G.nodes
                outfile.write('if_precursor_in_biggest_island,{}\n'.format(ifPrecursorInMain))
                outfile.write('if_target_and_precursor_in_the_same_island,{}\n'.format(self.ifTargetAndPrecursorInSameIsland()))
            outfile.write('number of reactions total,{}\n'.format(self.countNumReactions()))
        outfile.close()


    def countNumReactions(self):
        reactions_total = []
        for (u, v, reactions) in self.G.edges.data('reactions'):
            reactions_total.extend(reactions)
        return len(set(reactions_total))

    def ifTargetAndPrecursorInSameIsland(self):
        for island in list(self.G.subgraph(c) for c in nx.connected_components(self.G)):
            if self.data.main_precursor in island.nodes() and self.data.main_target in island.nodes():
                return True
        return False

    def getAllPathwaysRemainder(self):
        allPathways = nx.all_simple_paths(self.G, self.data.main_precursor, self.data.main_target)
        print(len(list(allPathways)))

    # ===============================================================================
    #
    # Draw graph for visual inspection ---------->

    def drawSubnetworkGraph(self, name):
        if len(self.G) > 0:
            pos=nx.spring_layout(self.G)
            nx.draw_networkx_edges(self.G, pos = pos)
            nx.draw_networkx_nodes(self.G, pos = pos, node_size=2, node_color='silver', node_shape='o')
            model_nodes_in_graph = list(set(self.G.nodes).intersection(self.data.modelMetabolites))
            nx.draw_networkx_nodes(self.G, pos, nodelist=model_nodes_in_graph, node_size=5, node_color='steelblue', node_shape='o', alpha=None)
            nx.draw_networkx_nodes(self.G, pos, nodelist=[self.data.main_target], node_size=50, node_color='gold', node_shape='o', alpha=None)
            plt.savefig(self.data.fig_dir+'/'+name+".png")
            plt.clf()

    # ===============================================================================
    #
    # Save graph in simple table format for Gephi (for more refined pictures)

    def dumpGraphForGephi(self):
        graph_file = open(self.data.fig_dir+'/graph_file_for_gephi.tsv', 'w')
        for edge in self.G.edges():
            graph_file.write('{}\t{}\n'.format(edge[0], edge[1]))
        graph_file.close()

        outfile = open(self.data.fig_dir+'/graph_file_for_gephi.gdf', 'w')

        #processing nodes of the subgraph
        outfile.write("nodedef>name VARCHAR,label VARCHAR,labelvisible BOOLEAN,width DOUBLE,height DOUBLE,color VARCHAR\n")

        df_nodes = pd.DataFrame(list(self.G.nodes()), columns=['cUID'])
        df_nodes.astype({'cUID':str})
        cmps_df = self.data.df_compounds.merge(df_nodes, on = 'cUID', how='inner')

        try:
            for node in self.G.nodes():
                cmp_df = cmps_df[cmps_df['cUID']==node].iloc[0]
                if node == self.data.main_target:
                    labvisib = 'true'
                    size = 9
                    color = "255, 215, 0"
                elif node in self.data.modelMetabolites:
                    labvisib = 'true'
                    size = 7
                    color = "0, 100, 0"
                else:
                    labvisib = 'false'
                    size = 5
                    color = "128, 128, 1"
                name=str(cmp_df['COMMON_NAME']).replace("'","")


                outfile.write("{},'{}','{}',{},{},'{}'\n".format(node, name, labvisib, size, size, color))

            #processing edges of the subgraph
            outfile.write("edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN, color VARCHAR\n")

            for edge in self.G.edges():
                outfile.write("{},{},{},'{}'\n".format(edge[0], edge[1], 'false','47, 79, 79'))

            outfile.close()

        except Exception as e:
            print('There was as exeption when producting Gephi file')
            print(e)


