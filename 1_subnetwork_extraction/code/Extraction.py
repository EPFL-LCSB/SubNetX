__author__ = 'anastasia'

# ===============================================================================
#
#  Expansion all
#
# ===============================================================================
import os
import time
from datetime import datetime
import multiprocessing

from Subnetwork import *

class Extraction():
    def __init__(self, data, graph):
        """

        :param data: data of the whole run: dataframes, network, parameters, lists
        :return:
        """
        self.data = data
        self.graph = graph


    def run(self):
        """
        Execution of the extraction called from Main.py
        :return:
        """
        # file to collect stats about how many pathways found at each round, how many boundary, what were the parameters etc
        self.openStatsFile()

        # define subnetwork which is going to be filled with linear pathways found towards precursor and towards the model
        self.subnetwork = Subnetwork(self.data)

        # number of initial pathways (shortest pathways) determined in the parameters
        self.graph.setNumberOfPathways(self.data.num_shortest_pathways)
        self.graph.setDistanceType(self.data.distance_transformation)

        # find initial pathways set
        # length of the longest pathway of the set determines the expansion limit for the expansion of the subnetwork
        pairsHaveRxns = self.getIntialPathwaysSet() # returns True if all pairs were successfully annotated with reactions

        if not pairsHaveRxns:
            print('Error: could not run algorithm since not all pairs are assigned reactions! Initial pathway search stage')

        #if no expansion the script can be used to just extract pathways to the model
        if self.data.run_expansion == True and pairsHaveRxns:
            # set number of pathways found towards each model precursor for all boundary metabolites
            self.graph.setNumberOfPathways(self.data.num_pathways_to_model)

            ##########    Expansion starts   ############
            flag_done = False
            while flag_done == False:
                self.data.current_round_extraction+=1
                flag_done = self.expansionRoundMain()

            ###########    Expansion ends    ############
            self.closeStatsFile()

        else: self.stats_file.close()

        self.subnetwork.writeGraphGpickle(self.data.auxilary_output_dir+'/subnetworkExtracted.gpickle')

    # ===============================================================================
    #
    #  Expansion functions

    # ===============================================================================
    #
    #  Generate initial pathways set

    def getIntialPathwaysSet(self):

        # if precursor defined, search towards the precursor. If no, search towards the model
        if self.data.main_precursor != 'all':
            initPathways = self.findInitialPathwaysSet() # and add to the subnetwork

        else:
            initPathways = self.findPathwaysForAllBoundary([self.data.main_target])
            # target is treated as the only boundary in round 0

        self.appendSubnetwork(initPathways)

        if self.data.run_expansion == True:
            return self.annotatePairs() # returns True if all pairs wre successfully annotated with reactions


    def findInitialPathwaysSet(self):
        t2 = time.time()
        initialPathways = self.graph.findInitialPathways(
            self.data.main_precursor, self.data.main_target, self.data.num_shortest_pathways)
        t3 = time.time()

        self.stats_file.write('{},'.format(len(initialPathways))) # number of pathways found for stats
        self.stats_file.write('{},'.format(t3 - t2)) # time to find initial pathways set

        # write down the initial pathways set to check later if it is recovered in the final network
        init_pathways_file = open(self.data.NICEpathways, 'w')
        for path in initialPathways:
            init_pathways_file.write('{}\n'.format('->'.join([str(i) for i in path])))
        init_pathways_file.close()

        return initialPathways

    # ===============================================================================
    #
    # Graph search towards the selected model ---------->

    def findPathwaysForAllBoundary(self, total_boundary):
        """
        Execute linear pathway search for all boundary that are not in the graph or the model
        :param total_boundary:
        :return:
        """

        ts0 = time.time()
        PathwaysList = []

        # pathway search is parallelized
        pool = multiprocessing.Pool()
        PathwaysList.extend(pool.map(self.graph.findShortestPathwaysToModel, list(total_boundary)))
        # pool has to be closed to avoid OSError: [Errno 24] Too many open files
        pool.close()

        # without multiprocess -->
        #for boundary in total_boundary:
        #    print('Boundary', boundary)
        #    PathwaysList.extend(self.graph.findShortestPathwaysToModel(boundary)) # <--

        ts1 = time.time()
        print('Time to find all pathways: ', ts1 - ts0)
        totalPathwaysList = []
        for path in PathwaysList:
            totalPathwaysList.extend(path)

        self.stats_file.write('{},'.format(len(totalPathwaysList))) # number of pathways found for stats
        self.stats_file.write('{},'.format(ts1 - ts0)) #  time to find all pathways for stats

        # write down the initial pathways set to check later if it is recovered in the final network
        # if the file does not exist means here we found an initial pathways set without defined precursor
        if not os.path.exists(self.data.NICEpathways):
            init_pathways_file = open(self.data.NICEpathways, 'w')
            for path in totalPathwaysList:
                init_pathways_file.write('{}\n'.format('->'.join([str(i) for i in path])))
            init_pathways_file.close()

        return totalPathwaysList


    def appendSubnetwork(self, totalPathwaysList):
        """
        adding all pathways to subnetwork
        :param totalPathwaysList: list of pathways found with linear pathway search
        :return: adds the pathway to the subnetwork
        """
        for path in totalPathwaysList:
            self.subnetwork.addPathwayFromPathwaySearchOutputToSubnetwork(path)

    def annotatePairs(self):
        """
        annotating pairs with reactions and boundaries
        :return:
        """
        ts2 = time.time()
        self.subnetwork.getPairsToAnnotate()
        self.stats_file.write('{},'.format(len(self.subnetwork.pairs_to_annotate))) # number of pairs to annotate for stats

        if not self.subnetwork.annotateNewPairsWithReactionsAndBoundary():
            return False # if some pair does not have any reaction
        ts3 = time.time()
        self.stats_file.write('{},'.format(ts3 - ts2)) # time to annotate pairs for stats
        return True

    # ===============================================================================
    #
    #  Expansion round

    def expansionRoundMain(self):
        """
        Round:
        1. Identify boundaries
        2. Find pathways toward the model
        :return:
        """

        # Identifying boundary

        # all total_boundary
        total_boundary_round = self.subnetwork.identifyOuterBoundaries()
        self.stats_file.write('{},'.format(len(total_boundary_round)))


        # total_boundary for which we can find pathways within the graph
        # IMPORTNANT: identify boundaries that are searcheable, meaning they are in the graph and not in the
        # list of boundaries for which the pathway was previously not possible to find
        total_boundary_searchable_round = set(total_boundary_round).intersection(set(self.graph.G.nodes))\
                                          -set(self.data.boundary_no_path)-set(self.data.boundaries_total)

        # Boundaries for which pathway cannot be found since the boundary is not in the graph
        total_boundary_unsearchable_round = set(total_boundary_round)-set(self.graph.G.nodes)\
                                          -set(self.data.boundary_no_path)-set(self.data.boundaries_total)
        self.data.boundaries_total.extend(total_boundary_searchable_round)
        self.stats_file.write('{}\n'.format(len(total_boundary_searchable_round)))

        # write down the boundaries of this round of extraction
        with open(self.data.boundaries_dir+'/round'+str(self.data.current_round_extraction)+'searchable.txt', 'w') as w:
            for boundary in total_boundary_searchable_round:
                w.write('{}\n'.format(boundary))

        with open(self.data.boundaries_dir+'/round'+str(self.data.current_round_extraction)+'unsearchable.txt', 'w') as w:
            for boundary in total_boundary_unsearchable_round:
                w.write('{}\n'.format(boundary))

        # main part - finding the pathways towards the model
        print('Boundaries in the graph:', len(total_boundary_searchable_round))
        if len(total_boundary_searchable_round) > 0:
            flag_done = False
            totalPathways = self.findPathwaysForAllBoundary(total_boundary_searchable_round)
            self.appendSubnetwork(totalPathways)
            pairsHaveRxns = self.annotatePairs() # returns True if all pairs wre successfully annotated with reactions
            if not pairsHaveRxns:
                exit('Error: could not run algorithm since not all pairs are assigned reactions! Expansion stage')

            with open(self.data.boundaries_dir+'/round'+str(self.data.current_round_extraction)+'_pathways.txt', 'w') as w:
                for pathway in totalPathways:
                    w.write('{}\n'.format('->'.join(pathway)))

        else:
            flag_done = True
            print("Pathways to model found for all boundary")

        return flag_done

    # ===============================================================================
    #
    # Opening the stats file, avoiding to overwrite the previous stats ---------->

    def openStatsFile(self):
        if not os.path.exists(self.data.stats_dir+'/subnetwork_extraction_stats.csv'):
            self.stats_file = open(self.data.stats_dir+'/subnetwork_extraction_stats.csv', 'w')
            self.stats_file.write('Pathways found,'\
    'Pathway search time,'\
    'Pairs to annotate,'\
    'DB annotation time,'\
    'Boundaries found total,'\
    'Boundaries found in the graph\n')
        else:
            self.stats_file = open(self.data.stats_dir+'/subnetwork_extraction_stats.csv', 'w')
            self.stats_file.write('iterating on {}\n'.format(datetime.now().isoformat(timespec='minutes')))

    def closeStatsFile(self):

        # write the final stats for the extracted network
        total_boundary = self.subnetwork.identifyOuterBoundaries()
        self.stats_file.write('{},'.format(len(total_boundary)))

        total_boundary_searchable = set(total_boundary).intersection(set(self.graph.G.nodes))

        self.stats_file.write('{}\n'.format(len(total_boundary_searchable)))

        #close the file
        self.stats_file.close()
