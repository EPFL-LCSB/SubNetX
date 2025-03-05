__author__ = 'anastasia'
##
## @file   Main.py
## @brief  Extracts a subnetwork based on the linear pathway search, reactions and boundary metabolites
##
## <!--------------------------------------------------------------------------
## This program follows the following logic:
## 1. Finds set of linear shortest pathways from the predefined precursor to predefined target
## 2. Assigns reactions to each step of the pathway, identifies boundary metabolites
## 3. Extracts the subnetwork around the initial set of pathways towards the model from ATLAS graph until expansion
##    limit is reached
## 4. Converges the subnetwork: removes edges that require metabolites outside the model and outside the subnetwork
##    until there are no more edges to remove
## 5. Evaluates in the precursor and target are located in the same graph component
##
## ------------------------------------------------------------------------ -->
##
import os

import sys
# Get the script's directory and add it to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

from Data import Data
from Extraction import *
from Convergence import *
from Compound import *
from Balance import *
from Format import *
from Subnetwork import *

def main():

    # get the data: pathway search parameters, list of organism metabolites, rxn and compounds dataframe
    data = Data()

    ####################################################################################################################
    #              BALANCE : crucial for correct optimization and yield calculations                                   #
    ####################################################################################################################
    bal = Balance(data)
    if not os.path.exists(data.reaction_balance):
        bal.createBalanceFile()

    ####################################################################################################################
    #              STRUCTURE BASED COMPOUNDS FILTERING :                                                               #
    #              1. helps to avoid nonsensical routes                                                                #
    #              2. helps to minimize the pathway search time                                                        #
    ####################################################################################################################
    cf = compoundFiltering(data) # initiating class for compound filtering
    data.precursorCompounds = set(cf.getCompoundsWithCarbon(data.modelMetabolites)) # filtering out the precursors without carbon
    print('Model metabolites with Carbon:', len(data.precursorCompounds))

    # create empty graph object
    graph = Graph(data)

    # extract total network from ATLAS files and transform it according to the parameters
    graph.formG()

    # identify set of model metabolites which can serve as precursors in graph search
    # setting data to have the model info the graph
    # Only the lists of metabolites in ecoli and yeast are used for pathway search
    graph.findPrecursorCompoundsInGraph()
    data.precursorCompounds = graph.precursorCompounds

    filterP = filterPrecursors(data) # init class from Compound file

    if data.filter_precursor_structure == True:
        filtered_metabolites_list = filterP.getListMatchingCompounds(data.compound_name, data.precursorCompounds)
        data.precursorCompounds = set(filtered_metabolites_list)

    ####################################################################################################################
    #              RUN EXPANSION :                                                                                     #
    #              1. find routes                                                                                      #
    #              2. find boundaries                                                                                  #
    #              3. find routes to boundaries                                                                        #
    #              4. repeat until the end of reaction network in reached                                              #
    ####################################################################################################################


    # run extraction of the subnetwork from the ATLAS graph and calculate the statistics on rounds
    e = Extraction(data, graph)
    e.run()

    # get the subnetwork that converged after number of rounds and calculate the statistics
    c = Convergence(data)
    c.run()

    ####################################################################################################################
    #                                             BOUNDARIES                                                           #
    ####################################################################################################################
    # write all boundaries by generation for the later analysis
    data.writeOutputBoundaryInfoList()

    ####################################################################################################################
    #                                    PRODUCE INPUT FOR OPTIMIZATION STEP                                           #
    ####################################################################################################################
    # write reactions LCSBID of the converged subnetwork for the next step: conversion to optimization
    if len(c.subnetwork.G) > 0:
        data.writeReactionsAll(c.subnetwork.G)

        f = Format(data)
        f.produceFormattedOutput()
    else:
        print("No balanced pathways possible!")

if __name__ == '__main__':
     main()



