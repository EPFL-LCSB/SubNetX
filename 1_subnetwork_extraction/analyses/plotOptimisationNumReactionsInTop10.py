__author__ = 'anastasia'

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import os

# in case needed to do for several modes at the same time
#folders = [i for i in os.listdir('data') if not i.startswith('.')]

folders = ['biosubnet_results_expanded_mode']

compounds = [i for i in os.listdir('data/'+folders[0]) if not (i.startswith('.') or i=='scopolamine')]

def numReactions(row):
    return len(row['reactions'].split(' // '))

def plotLine():
    """ First attempt to plot - Vassily did not like this representaiton """
    legend_lines = []
    legend_compounds = []
    ax = plt.figure().gca()
    for cmp in compounds:
        if not os.path.exists('data/'+folders[0]+'/'+cmp+'/pthws_ids_1_min.csv'):
            continue
        df=pd.read_csv('data/'+folders[0]+'/'+cmp+'/pthws_ids_1_min.csv')
        #print(df.columns) #['Pathways found', 'Pathway search time', 'Pairs to annotate', 'DB annotation time', 'Boundaries found total','Boundaries found in the graph']
        df['round'] = range(1, len(df) + 1)
        df['numReactiions'] = df.apply(numReactions, axis=1)
        print(df['numReactiions'].to_list())
        line, = ax.plot(df['round'].to_list(), df['numReactiions'].to_list(),'-o')
        legend_lines.append(line)
        legend_compounds.append(cmp)

    ax.legend(legend_lines, legend_compounds)
    ax.set_xlabel('Rank of the pathway')
    ax.set_ylabel('Number of reactions')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig('figures/optimization_reactions_num_'+folders[0]+'.eps', format='eps', bbox_inches='tight')
    plt.close()

def getBars(listValues):
    values = list(set(listValues))
    counts = []
    for v in values:
        counts.append(listValues.count(v))
    return values, counts

def plotBar():
    colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
    num_axes = 0
    for index, cmp in enumerate(compounds):
        if os.path.exists('data/'+folders[0]+'/'+cmp+'/pthws_ids_1_min.csv'):
            num_axes+=1

    #fig, axs = plt.subplots(1, num_axes)
    ax = plt.figure(figsize=(30,10)).gca()
    for index, cmp in enumerate(compounds):
        legend_lines = []
        legend_compounds = []

        if not os.path.exists('data/'+folders[0]+'/'+cmp+'/pthws_ids_1_min.csv'):
            continue

        #ax = axs[index]
        df=pd.read_csv('data/'+folders[0]+'/'+cmp+'/pthws_ids_1_min.csv')
        df['numReactiions'] = df.apply(numReactions, axis=1)
        values, counts = getBars(df['numReactiions'].to_list())
        line = ax.bar([0.4*index-0.8 + v if index>1 else v for v in values ], counts, 0.4, color = colors[index])
        legend_lines.append(line)
        legend_compounds.append(cmp)

    ax.legend(legend_lines, legend_compounds)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax.set_xlabel('Number of reactions in a pathway', size=30)
    ax.set_ylabel('Number of pathways', size=30)
    #ax.set_ylim([0, 12])

    ax.set_xticks([5, 6, 7, 21, 22, 24])

    plt.xticks(size = 30)
    plt.yticks(size = 30)
    #plt.show()

    plt.savefig('figures/optimization_reactions_bar_num_'+folders[0]+'.eps', format='eps', bbox_inches='tight')
    plt.close()

plotBar()