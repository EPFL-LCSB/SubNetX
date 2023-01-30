__author__ = 'anastasia'
from Data import Data
from Extraction import *
from Convergence import *

class Tests():

    # Generate subnetwork for each individual pathway and check convergence --------------------------------------- -->
    def testIndividualPathwaysSubnetworks(self):

        # get the data: pathway search parameters, list of organism metabolites, rxn and compounds dataframe
        self.data = Data()
        self.data.readParametersFile()
        self.data.getModelMetabolites()
        self.data.getReactionsAndBoundaryDataframe()
        self.data.createDirectories() # create directories for the output

        self.extr = Extraction(self.data)
        self.conv = Convergence(self.data)

        # file to collect stats about how many pathways found at each round, how many boundary, what were the parameters etc
        self.extr.openStatsFile()

        # find set of initial pathways
        initPathways = self.extr.findInitialPathwaysSet()

        # set number of pathways found towards each model precursor for all boundary metabolites
        self.extr.graph.setNumberOfPathways(self.extr.data.num_pathways_to_model)

        # for each pathway extract and converge the subnetwork
        #initPathways = [initPathways[0], initPathways[4], initPathways[-1]]
        initPathways = [[1467915486, 1142614524, 1142614546, 1467892230, 1467886235, 1468725650, 1467900788, 1467883653, 1467874237]]
        for pwRank, pw in enumerate(initPathways):

            # define subnetwork which is going to be filled with linear pathways found towards precursor and towards the model
            self.extr.subnetwork = Subnetwork(self.data)

            # add only this pw to the initial subnetwork
            print('Current PW:', pwRank, pw)

            self.extr.subnetwork.addPathwayFromPathwaySearchOutputToSubnetwork(pw)

            #self.extr.subnetwork.getPairsToAnnotate()
            #self.extr.subnetwork.annotateNewPairsWithReactionsAndBoundary(self.extr.data.dfCmp1Cmp2RxnCompounds)
            #total_boundary = self.extr.subnetwork.identifyOuterBoundaries(self.data.modelMetabolites)
            #print(total_boundary)
            #del self.extr.subnetwork
            #continue

            print('Start extraction')
            self.extractSubnetworkPerPathway()

            # Reassign subnetwork from extraction to convergence
            self.conv.subnetwork = self.extr.subnetwork

            print('Test convergeance of pathway')
            self.testConvergeancePerPathway(pwRank)
            print()

            # close the stats file of network extraction if all the pathways are analysed
            if pwRank == len(initPathways)-1:
                self.extr.closeStatsFile()

            del self.extr.subnetwork
            del self.conv.subnetwork

    def extractSubnetworkPerPathway(self):

        ts2 = time.time()
        self.extr.subnetwork.getPairsToAnnotate()
        self.extr.stats_file.write('{}\t'.format(len(self.extr.subnetwork.pairs_to_annotate))) # number of pairs to annotate for stats

        self.extr.subnetwork.annotateNewPairsWithReactionsAndBoundary(self.extr.data.dfCmp1Cmp2RxnCompounds)
        ts3 = time.time()
        self.extr.stats_file.write('{}\t'.format(ts3 - ts2)) # time to annotate pairs for stats

        #  Expansion
        flag_done = False
        count_round = 0
        while flag_done == False:
            count_round += 1
            print('Current round:', count_round)
            flag_done = self.extr.expansionRoundMain()

    def testConvergeancePerPathway(self, pwRank):

        self.conv.subnetwork.calculateStatistics('converged_graph_extracted_stats_pw{}.txt'.format(pwRank))
        self.conv.openStatsFileConvergeance('covergeance_rounds_stats_pw{}.tsv'.format(pwRank))

        # Check if pathways converged after number of rounds ---------->
        ifConverged = False
        while ifConverged == False:
            ifConverged = self.conv.convergeNetwork()

        self.conv.stats_file_convergeance.close()
        self.conv.subnetwork.calculateStatistics('converged_graph_stats_pw{}.txt'.format(pwRank))

    def howManyInitialPathwaysInSubnetwork(self):
        data = Data()
        data.readParametersFile()
        G = nx.read_gpickle('output/'+data.compound_name+'/subnetworkConverged.gpickle')
        init_pathways = data.getInitialPathways()
        pw_converged = []
        pw_not_converged = []
        for pw in init_pathways:
            flag = True
            for i in range(len(pw)-1):
                u = pw[i]
                v = pw[i+1]
                if not ((u, v) in G.edges() or (v, u) in G.edges()):
                    flag = False
            if flag == False:
                pw_not_converged.append(pw)
            else:
                pw_converged.append(pw)
        print('Pathways that did not converge: ', len(pw_not_converged))
        with open('output/'+data.compound_name+'/initial_pathways_converged.txt', 'w') as w:
            for pw in pw_converged:
                w.write('{}\n'.format('->'.join(str(i) for i in pw)))



tests = Tests()
tests.testIndividualPathwaysSubnetworks()
#tests.howManyInitialPathwaysInSubnetwork()

# Functions for drawing the subnetwork
def drawSubnetworkExtracted(data):
    data=Data()
    data.readParametersFile()
    data.getModelMetabolites()
    subnetwork = Subnetwork(data)
    G = subnetwork.getGraph(data.auxilary_output_dir+'/subnetworkExtracted.gpickle')
    subnetwork.G = G
    subnetwork.drawSubnetworkGraph("extractedSubnetwork")

def dumpSubnetworkExtracted():
    data=Data()
    data.readParametersFile()
    data.getModelMetabolites()
    subnetwork = Subnetwork(data)
    G = subnetwork.getGraph(data.auxilary_output_dir+'/subnetworkExtracted.gpickle')
    subnetwork.G = G
    subnetwork.dumpGraphForGephi()

#drawSubnetworkExtracted()
#dumpSubnetworkExtracted()