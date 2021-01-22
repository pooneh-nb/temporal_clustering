import operator

import networkx as nx
import json
from statistics import mean
import matplotlib.pyplot as plt

import pandas as pd

def create_graph(path):
    G = nx.Graph()
    for line in open(path):
        if line.find('#') and line.find('%') < 0 and line != '\n':
            line = line.strip('\n').split('\t')
            if int(line[0]) == int(line[1]):
                continue

            if G.has_edge(int(line[0]), int(line[1])):
                if int(line[2]) not in G.edges[int(line[0]), int(line[1])]['t']:
                    G.edges[int(line[0]), int(line[1])]['t'].append(int(line[2]))
                    G.edges[int(line[0]), int(line[1])]['w'].append(int(line[3]))
                    # print(int(line[0]), int(line[1]) ,G.edges[int(line[0]), int(line[1])])
            else:
                G.add_edge(int(line[0]), int(line[1]),  t=[int(line[2])], w=[int(line[3])])
    return G

def masure_avg_weight():
    path = "temporal_df_2016"
    FG = create_graph(path)
    #FG = nx.Graph()
    #FG.add_weighted_edges_from([(1, 2, 0.125), (1, 3, 0.75), (2, 4, 1.2), (3, 4, 0.375)])
    degree_sequence = sorted([d for n, d in FG.degree()], reverse=True)
    dmax = max(degree_sequence)
    print(dmax)

    wt = {}
    for n, nbrs in FG.adj.items():
        wt[n] = 0
        for nbr, eattr in nbrs.items():
            wt[n] = wt[n] + round(mean(eattr['w']), 2)

    weight_max = max(wt.items(), key=operator.itemgetter(1))[0]
    print(wt[weight_max])

    #with open('weights.json', 'wt') as fp:
        #json.dump(wt, fp, indent=4)


def weight_histogram():
    with open('weights.json', 'rt') as ww:
        weights = json.load(ww)
    weight_data = []
    for key, weight in weights.items():
        if 5000<weight < 100000:
         weight_data.append(round(weight,2))
    print(len(weight_data))

    fig, axs = plt.subplots(1, 1)
    axs.hist(weight_data, bins=10)
    plt.show()



#masure_avg_weight()
#weight_histogram()