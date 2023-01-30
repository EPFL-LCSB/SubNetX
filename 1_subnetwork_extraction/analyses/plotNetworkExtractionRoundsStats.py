__author__ = 'anastasia'

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import os

def plotExtraction(fold):
    compounds = [i for i in os.listdir('data/'+folders[0]) if not (i.startswith('.'))]
    legend_lines = []
    legend_compounds = []
    ax = plt.figure().gca()
    for cmp in compounds:
        df=pd.read_csv('data/'+fold+'/'+cmp+'/stats/subnetwork_extraction_stats.csv')
        #print(df.columns) #['Pathways found', 'Pathway search time', 'Pairs to annotate', 'DB annotation time', 'Boundaries found total','Boundaries found in the graph']
        df['round'] = range(1, len(df) + 1)
        line, = ax.plot(df['round'].to_list(), df['Boundaries found total'].to_list(),'-o')
        legend_lines.append(line)
        legend_compounds.append(cmp)
    ax.legend(legend_lines, legend_compounds)
    ax.set_xlabel('Round of extraction')
    ax.set_ylabel('Number of boundaries')
    plt.savefig('figures/subnetwork_extraction_rounds_'+fold+'.eps', format='eps', bbox_inches='tight')
    plt.close()

def plotConvergeance(fold):
    compounds = [i for i in os.listdir('data/'+folders[0]) if not (i.startswith('.'))]
    legend_lines = []
    legend_compounds = []
    ax = plt.figure().gca()
    for cmp in compounds:
        df=pd.read_csv('data/'+fold+'/'+cmp+'/stats/covergeance_rounds_stats.csv')
        df['round'] = range(1, len(df) + 1)
        line, = ax.plot(df['round'].to_list(), df['Number of all outer boundaries'].to_list(),'-o')
        legend_lines.append(line)
        legend_compounds.append(cmp)
    ax.legend(legend_lines, legend_compounds)
    ax.set_xlabel('Round of extraction')
    ax.set_ylabel('Number of boundaries')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig('figures/subnetwork_convergeace_rounds_'+fold+'.eps', format='eps', bbox_inches='tight')
    plt.close()

def plotTotalConvergedStats(fold):
    compounds = [i for i in os.listdir('data/'+fold) if not (i.startswith('.'))]
    legend_lines = []
    legend_compounds = []
    #fig12, (ax1, ax2) = plt.subplots(2, 1)
    #fig34, (ax3, ax4) = plt.subplots(2, 1)

    listDF = list()
    for cmp in compounds:
        statsfile = [i for i in os.listdir('data/'+fold+'/'+cmp+'/stats/') if i.startswith('converged_graph_stats')][0]
        dict_stats = dict()
        with open('data/'+fold+'/'+cmp+'/stats/'+statsfile) as f:
            for line in f:
                dict_stats['compound_name']=cmp
                dict_stats[line.split(',')[0]]=line.strip().split(',')[1]
            listDF.append(dict_stats)
    df = pd.DataFrame(listDF)

    #line, = ax1.plot(df['compound_name'].to_list(), df['Number of all outer boundaries'].to_list(),'-o')
    # legend_lines.append(line)
    # legend_compounds.append(cmp)
    # ax.legend(legend_lines, legend_compounds)
    # ax.set_xlabel('Round of extraction')
    # ax.set_ylabel('Number of boundaries')
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # plt.savefig('figures/subnetwork_convergeace_rounds_'+fold+'.eps', format='eps', bbox_inches='tight')
    # plt.close()

    #df.to_csv('output_convergence_total_stats.csv')


#folders = [i for i in os.listdir('data') if not i.startswith('.')]
folders = ['biosubnet_results_expanded_mode']

for fold in folders:
    plotExtraction(fold)
    plotConvergeance(fold)
#plotTotalConvergedStats(folders[0])