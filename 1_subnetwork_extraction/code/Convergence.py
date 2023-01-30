__author__ = 'anastasia'
# ===============================================================================
#
# Convergence check for the extracted subnetwork
#
# ===============================================================================

from Subnetwork import *

# Class convergence removes all the connections that are not balanced

class Convergence():
    def __init__(self, data):
        self.data = data

    def run(self):
        """
        Function called from Main.py to execute convergence
        :return: produces subnetwokConverged.gpickle, stats of convergence
        """
        self.subnetwork = Subnetwork(self.data)

        self.subnetwork.G = self.subnetwork.getGraph(self.data.auxilary_output_dir+'/subnetworkExtracted.gpickle')

        self.subnetwork.calculateStatistics()
        self.openStatsFileConvergeance()

        # Check if pathways converged after number of rounds ---------->
        ifConverged = False
        while ifConverged == False:
            ifConverged = self.convergeNetwork()

        self.stats_file_convergeance.close()
        self.subnetwork.calculateStatistics()
        self.subnetwork.drawSubnetworkGraph("convergedIslandsDeleted_nodes{}".format(len(self.subnetwork.G.nodes())))
        self.subnetwork.writeGraphGpickle(self.data.auxilary_output_dir+'/subnetworkConverged.gpickle')
        self.subnetwork.dumpGraphForGephi()

    def convergeNetwork(self):
        """
        Rounds of convergence:
        remove non converged edges and write stats
        :return:
        """

        # remove edges that did not converge
        self.removeNonConvergedEdges()
        # Write stats
        self.stats_file_convergeance.write('{},'.format(len(self.subnetwork.G.edges())))
        self.stats_file_convergeance.write('{},'.format(len(self.subnetwork.G.nodes())))

        # delete islands that are not connected to the model
        ifConverged = self.deleteIslandsWithoutModelMetabolites()
        self.stats_file_convergeance.write('{}\n'.format(len(self.subnetwork.G.nodes())))

        return ifConverged

    def removeNonConvergedEdges(self):
        """
        Identify edges that have boundaries outside the graph and delete them from the graph
        :return:
        """

        # all outer boundary metabolites besides the model metabolites and metabolites in the graph
        allOuterBoundaries = self.subnetwork.identifyOuterBoundaries()
        # write stats
        self.stats_file_convergeance.write('{},'.format(len(allOuterBoundaries)))
        self.stats_file_convergeance.write('{},'.format(len(self.subnetwork.G.edges())))
        # collect edges to remove : needed since graph cannot be parsed and modified at the same time
        edges_remove = []
        for edge in self.subnetwork.G.edges():
            if any(i in allOuterBoundaries for i in self.subnetwork.G[edge[0]][edge[1]]['all_rxns_boundary']):
                edges_remove.append(edge)
        self.stats_file_convergeance.write('{},'.format(len(edges_remove)))
        # remove edges
        for edge in edges_remove:
            self.subnetwork.G.remove_edge(edge[0], edge[1])

    def deleteIslandsWithoutModelMetabolites(self):
        """
        Delete network components that are not connected to any model metabolites -> cannot be produced
        :return:
        """
        remove_nodes = []
        for island in list(self.subnetwork.G.subgraph(c) for c in nx.connected_components(self.subnetwork.G)):
            if not (any(i in island.nodes() for i in self.data.modelMetabolites) or (self.data.main_target in island.nodes())):
                remove_nodes.extend(island.nodes)
        for node in remove_nodes:
            self.subnetwork.G.remove_node(node)
        ifConverged = len(remove_nodes)<1
        return ifConverged

    def openStatsFileConvergeance(self):
        """
        Initiate the file to write the convergence statistics
        :return:
        """
        self.stats_file_convergeance = open(self.data.stats_dir+'/covergeance_rounds_stats.csv', 'w')
        self.stats_file_convergeance.write('Number of all outer boundaries,Edges in subnetwork before removing the edges with nonconverged boundary,Edges to remove,'\
    'Edges in subnetwork after removing the edges with nonconverged boundary,'\
    'Nodes in subnetwork before removing the nodes on islands without model metabolites,'\
    'Nodes in subnetwork after removing the nodes on islands without model metabolites\n')
