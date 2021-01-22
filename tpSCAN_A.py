# -*- coding: utf-8 -*-
import readfile
import math
import sys
import os
import matplotlib.pyplot as plt
import networkx as nx
from time import sleep, time
import datetime
import json
from statistics import mean
import math
import gc
import itertools
import pprint
import multiprocessing
from multiprocessing import Pool as ThreadPool
import itertools


class tGraph:

    def __init__(self, path):
        self.path = path
        # list of (node, degree) sorted based on degree
        self.rank = []
        # the highest rank of rank list: (92,181) based on node's degree. here node 92 with degree 181 has the
        # highest rank
        self.toprank = 0
        self.visited_set = []
        self.visited_node = set() # add nodes which their eligibility for being a core has been tested
        self.visited_edge = set()
        # degree of similarity between two neighbor nodes
        self.eps = 0.3
        # minimum number of similar neighbors
        self.miu = 5
        # minimum number of similar timestamps
        self.tau = 3
        self.theta = 1
        self.sigma = {}
        # calculated sigma for each edge at timestamp t
        self.sigma_t = {}
        self.union_set = []  # list of cores
        self.collections = {}
        # dictionary of each node as key along with their neighbors and degrees as value
        #     {u: {v:v_deg, h:h_degree,...}
        #      v: {list of v's neighbors along with their degree}
        #     }
        self.adj = {}
        self.subgraph = {}
        self.frquent_set = {}
        # create a temporal graph in this format (u,v,t)
        self.G = readfile.tGraph(self.path)

        print(len(self.G.nodes()))
        print(len(self.G.edges()))
        temporal_edge = 0
        for nodex in self.G.nodes():
            for nodey in self.G.adj[nodex]:
                temporal_edge += len(self.G.edges[nodex, nodey]['t'])
        print("temporal edges:" + str(temporal_edge / 2))

        # dict: {node:degree}
        ranktemp = {}
        for node_temp in self.G.nodes():
            # add two attributes to each node; l and u. l=0 and u = number of neighbors(degree)
            self.G.nodes[node_temp]['l'] = 0
            # attr u should be the average weight of u
            #self.G.nodes[node_temp]['u'] = len(self.G.adj[node_temp])
            wt = 0
            for nbr, eattr in self.G.adj[node_temp].items():
                wt = wt + sum(eattr['w'])
            #self.G.nodes[node_temp]['u'] = len(self.G.adj[node_temp])
            self.G.nodes[node_temp]['u'] = wt
            ranktemp[node_temp] = self.G.nodes[node_temp]['u']
        # then sort ranktemp dictionary based on degree of each node
        # the node with largest degree is number one
        self.rank = sorted(ranktemp.items(), key=lambda item: item[1], reverse=True)

        for node_temp in ranktemp:
            self.adj[node_temp] = {}
            for item in self.G.adj[node_temp]:
                if item == node_temp:
                    print(item)
                else:
                    # the adj member is created here
                    self.adj[node_temp][item] = ranktemp[item]

            # for each node, sort its neighbors based on its degree
            adjtemp = sorted(self.adj[node_temp].items(), key=lambda item: item[1], reverse=True)

            self.adj[node_temp] = []
            for item in adjtemp:
                # add neighbors of each node, based on highest degree
                self.adj[node_temp].append(item[0])

        ranktemp = []
        # set to the node with the highest rank
        self.toprank = self.rank[0]
        for i in self.rank:
            ranktemp.append(i[0])
        self.rank = ranktemp

        del (ranktemp)

    def tDistribution(self, tempG):
        timestamps = {}
        for item in tempG.edges.data():
            if item[2]['t'] in timestamps:
                timestamps[item[2]['t']] = timestamps[item[2]['t']] + 1
            else:
                timestamps[item[2]['t']] = 1

        min = -1
        max = -1
        x = []
        y = []
        for k, v in timestamps.items():
            if k > max:
                max = k
            if min > k or min == -1:
                min = k
            x.append(k)
            y.append(v)
        print(min, max)
        print(time.gmtime(min), time.gmtime(max))

    def check_SCANA_core(self, u):
        # this function check whether node u can be a core. A node u E V is called
        # a (miu, tau, epsilon)-stable core if there exist a set of neighbors
        # N~(u) E N(u) of u that satisfies the following conditions:
        # 1. |N~(u)| >= miu
        # 2. there are at least tau snapshots containing the star-shaped structure formed by
        #    N~(u) of u
        # 3. in each snapshot, sigma(u,v) >= epsilon for any v E N~(u)

        # check whether the number of node u's neighbors is larger than miu

        # condition 1
        """if u in self.visited_node:
            return False"""
        initial_weight = 0
        if self.G.nodes[u]['l'] < self.miu and self.G.nodes[u]['u'] >= self.miu:
            candidate_set = []
            for v in self.G.adj[u]:
                # check all the neighbors of u. If any of them, (u,v) edge, has been repeated
                # in more tha tau snapshots then add that edge to candidate_set
                #print(self.G.edges[u, v]['t'])
                if len(self.G.edges[u, v]['t']) >= self.tau:
                    candidate_set.append(v)
                    #initial_weight = initial_weight + round(mean(self.G.edges[u,v]['w']),2)
                    initial_weight = initial_weight + sum(self.G.edges[u, v]['w'])
            # condition 2

            if initial_weight < self.miu:
            #if len(candidate_set) < self.miu:
                return False

            # find what is sum of the weights of candidate neighbors at each timestamp;
            candidate_time = {}
            for v in candidate_set:
                idx = 0
                for time_item in self.G.edges[u, v]['t']:
                    if time_item in candidate_time:
                        candidate_time[time_item] += self.G.edges[u,v]['w'][idx]
                        #candidate_time[time_item] += 1
                    else:
                        candidate_time[time_item] = self.G.edges[u,v]['w'][idx]
                        #candidate_time[time_item] = 1
                    idx += 1

            # only work on timestamps that their aggregated weight is more than miu
            times_more_than_miu = []
            for key in candidate_time:
                if candidate_time[key] >= self.miu:
                    times_more_than_miu.append(key)
            if len(times_more_than_miu) < self.tau:
                return False
            # condition 3

            # then check the similarity of each candidate edge at times_more_than_miu.
            # If the similarity is more than eps, add the weight of that edge to miu_calculate. Also update
            # the 'l' attribute of node u and its neighbors.
            # At the end of this method, if the self.G.nodes[u]['u'] is more than mui, it is eligible to be a core.
            tau_calculate = 0
            for t in times_more_than_miu:
                if tau_calculate >= self.tau:
                    break
                # sum of weights of desired neighbors at each timestamp
                miu_calculate = 0
                # sum of weights of all candidate neighbors regardless of their similarity
                max_miu_calculate = candidate_time[t]
                for v in candidate_set:
                    #if tau_calculate >= self.tau:
                        #break
                    if v <= u:
                        edge_set = (u, v, t)
                    else:
                        edge_set = (v, u, t)
                    if edge_set not in self.visited_edge:
                        if edge_set not in self.sigma_t:
                            # set the similarity of two nodes at time t
                            self.sigma_t[edge_set] = self.compute_sigma_at_one_time(u, v, t)

                        if self.sigma_t[edge_set] >= self.eps:
                            idx_t_u = self.G.edges[u, v]['t'].index(t)
                            miu_calculate += self.G.edges[u,v]['w'][idx_t_u]
                            self.G.nodes[u]['l'] += self.G.edges[u,v]['w'][idx_t_u]
                            self.G.nodes[v]['l'] += self.G.edges[u,v]['w'][idx_t_u]
                            self.visited_edge.add(edge_set)
                            #miu_calculate += 1
                        else:
                            if t not in self.G.edges[u, v]['t']:
                                pass
                            else:
                                idx_t_u = self.G.edges[u, v]['t'].index(t)
                                max_miu_calculate -= self.G.edges[u,v]['w'][idx_t_u]
                            #max_miu_calculate -= 1
                        #if miu_calculate >= self.miu:
                            #tau_calculate += 1
                            #continue

                        #if max_miu_calculate < self.miu:
                            #continue
                    else:
                        if v <= u:
                            edge_set = (u, v, t)
                        else:
                            edge_set = (v, u, t)
                        if edge_set not in self.sigma_t:
                            # set the similarity of two nodes at time t
                            self.sigma_t[edge_set] = self.compute_sigma_at_one_time(u, v, t)

                        if self.sigma_t[edge_set] >= self.eps:
                            idx_t_u = self.G.edges[u, v]['t'].index(t)
                            miu_calculate += self.G.edges[u, v]['w'][idx_t_u]
                        else:
                            if t not in self.G.edges[u, v]['t']:
                                pass
                            else:
                                idx_t_u = self.G.edges[u, v]['t'].index(t)
                                max_miu_calculate -= self.G.edges[u, v]['w'][idx_t_u]
                            # max_miu_calculate -= 1

                    if miu_calculate >= self.miu:
                        tau_calculate += 1
                        break
                        #continue

                    if max_miu_calculate < self.miu:
                        break
                        #continue

            if miu_calculate >= self.miu:
             tau_calculate += 1
             #continue

            if tau_calculate < self.tau:
                #self.G.nodes[u]['l'] = 0
                self.visited_node.add(u)
                return False

            """for v in candidate_set:
                edge_set = (u, v)
                if v >= u:
                    edge_set = (u, v)
                else:
                    edge_set = (v, u)

                if edge_set not in self.sigma:
                    self.sigma[edge_set] = self.compute_sigma(u, v)
                    # if the similarity of edge_set is more than sigma in tau snapshots
                    if self.sigma[edge_set] >= self.tau:
                        for t in times_more_than_miu:
                            if t in self.G.edges[u, v]['t']:
                                idx = self.G.edges[u,v]['t'].index(t)
                                self.G.nodes[u]['l'] += round(self.G.edges[u, v]['w'][idx], 2)
                                self.G.nodes[v]['l'] += round(self.G.edges[u, v]['w'][idx], 2)
                        self.G.nodes[u]['l'] /= len(times_more_than_miu)
                        self.G.nodes[v]['l'] /= len(times_more_than_miu)
                        #self.G.nodes[u]['l'] += 1
                        #self.G.nodes[v]['l'] += 1
                    else:
                        for t in times_more_than_miu:
                            if t in self.G.edges[u, v]['t']:
                                idx = self.G.edges[u, v]['t'].index(t)
                                self.G.nodes[u]['l'] -= round(self.G.edges[u, v]['w'][idx], 2)
                                self.G.nodes[v]['l'] -= round(self.G.edges[u, v]['w'][idx], 2)
                        #self.G.nodes[u]['u'] -= 1
                        #self.G.nodes[v]['u'] -= 1

                    if self.G.nodes[u]['l'] >= self.miu or self.G.nodes[u]['u'] < self.miu:
                        break"""

        self.visited_node.add(u)
        return False

    def cluster_SCANA_core(self, u):
        candidate_set = []
        for v in self.G.adj[u]:
            if len(self.G.edges[u, v]['t']) >= self.tau:
                candidate_set.append(v)

        for v in candidate_set:
            #edge_set = (u, v)
            if v <= u:
                edge_set = (u, v)
            else:
                edge_set = (v, u)

            if edge_set in self.sigma:
                if self.G.nodes[v]['l'] >= self.miu and self.sigma[edge_set] >= self.tau:
                    self.union(u, v)
            else:
                if self.G.nodes[v]['u'] >= self.miu:
                    self.sigma[edge_set] = self.compute_sigma(u, v)
                    """if self.sigma[edge_set] >= self.tau:
                        #self.G.nodes[u]['l'] += sum(self.G.edges[u, v]['w'])
                        self.G.nodes[v]['l'] += sum(self.G.edges[u, v]['w'])
                        #self.G.nodes[u]['l'] += round(mean(self.G.edges[u, v]['w']))
                        #self.G.nodes[v]['l'] += round(mean(self.G.edges[u, v]['w']))
                        #self.G.nodes[u]['l'] += 1
                        #self.G.nodes[v]['l'] += 1
                    else:
                        #self.G.nodes[u]['l'] -= sum(self.G.edges[u, v]['w'])
                        self.G.nodes[v]['l'] -= sum(self.G.edges[u, v]['w'])
                        #self.G.nodes[u]['u'] -= 1
                        #self.G.nodes[v]['u'] -= 1"""

                    if self.sigma[edge_set] >= self.tau:
                        self.check_SCANA_core(v)
                        if self.G.nodes[v]['l'] >= self.miu:
                            self.union(u, v)
                            #self.visited_node.add(u)

    def add_node_set(self, u):
        if len(self.union_set):
            for set in self.union_set:
                if u in set:
                    return 0
                if u not in set:
                    pass
            self.union_set.append([u])
            print(self.union_set)
        else:
            self.union_set.append([u])
            print(self.union_set)

    def union(self, u, v):
        if len(self.union_set):
            flag = 0
            set1 = []
            set2 = []
            for set in self.union_set:
                if u in set and v in set:
                    flag = -1  # no need to change
                    break
                if u in set and v not in set:
                    set1 = set
                    flag = flag + 1
                if v in set and u not in set:
                    set2 = set
                    flag = flag + 1
            if flag == 0:
                temp = [u, v]
                self.union_set.append(temp)
            if flag == 1:
                if set1:
                    index_temp = self.union_set.index(set1)
                    self.union_set[index_temp].append(v)
                if set2:
                    index_temp = self.union_set.index(set2)
                    self.union_set[index_temp].append(u)
            if flag == 2:
                self.union_set.remove(set1)
                self.union_set.remove(set2)
                union_temp = set1 + set2
                self.union_set.append(union_temp)
            if flag > 2:
                print("unnion error")

        else:
            temp = [u, v]
            self.union_set.append(temp)

    def compute_sigma_at_one_time(self, u, v, t):
        if t not in self.G.edges[u, v]['t']:
            return 0
        # ~~~~~~~~~~cosine similarity~~~~~~~~~~~~~~~
        w_adju = 0
        adju = []
        #idx_t_u = None
        for vertex in self.adj[u]:
            if t in self.G.edges[u, vertex]['t']:
                adju.append(vertex)
                idx_t_u = self.G.edges[u, vertex]['t'].index(t)
                w_adju += pow(self.G.edges[u, vertex]['w'][idx_t_u],2)

        w_adjv = 0
        adjv = []
        #idx_t_v = None
        for vertex in self.adj[v]:
            if t in self.G.edges[v, vertex]['t']:
                adjv.append(vertex)
                idx_t_v = self.G.edges[v, vertex]['t'].index(t)
                w_adjv += pow(self.G.edges[v, vertex]['w'][idx_t_v], 2)

        lenuadj = len(adju) + 1
        lenvadj = len(adjv) + 1
        if lenuadj < self.eps * self.eps * lenvadj or lenvadj < self.eps * self.eps * lenuadj:
            return 0

        common_nbr = set(adju).intersection(set(adjv))
        similarity_w = 0
        if len(common_nbr) == 0 and t in self.G.edges[u, v]['t']:
            idx_t_u = self.G.edges[u, v]['t'].index(t)
            similarity_w = pow(self.G.edges[u, v]['w'][idx_t_u], 2)
        else:
            idx_t_ = self.G.edges[u, v]['t'].index(t)
            similarity_w = pow(self.G.edges[u, v]['w'][idx_t_], 2)
            for nbr in common_nbr:
                idx_t_u = self.G.edges[u, nbr]['t'].index(t)
                idx_t_v = self.G.edges[v, nbr]['t'].index(t)
                similarity_w += self.G.edges[u, nbr]['w'][idx_t_u] * self.G.edges[v, nbr]['w'][idx_t_v]

        sigma = similarity_w / math.sqrt(w_adju*w_adjv)

        if sigma < self.eps:
            return 0
        else:
            return self.eps + 0.1
        """
        # list of neighbors of u at time t
        adju = []
        for vertex in self.adj[u]:
            if t in self.G.edges[u, vertex]['t']:
                adju.append(vertex)

        # list of neighbors of v at time t
        adjv = []
        for vertex in self.adj[v]:
            if t in self.G.edges[v, vertex]['t']:
                adjv.append(vertex)


        #lenuadj = len(adju) + 1
        lenuadj = w_adju
        #lenvadj = len(adjv) + 1
        lenvadj = w_adjv
        if lenuadj < self.eps * self.eps * lenvadj or lenvadj < self.eps * self.eps * lenuadj:
            # print(u,v,lenuadj,lenuadj)
            return 0

        len_v_u = len(set(adju) & set(adjv)) + 2
        if len_v_u < self.eps * math.sqrt(lenuadj * lenvadj):
            # print(u, v, len_v_u, lenuadj, lenuadj)
            return 0
        else:
            return self.eps + 0.1"""

    def compute_sigma(self, u, v):
        tau = 0

        if len(self.G.edges[u, v]['t']) < self.tau:
            return 0

        for t in self.G.edges[u, v]['t']:
            #edge_set = (u, v, t)
            if v <= u:
                edge_set = (u, v, t)
            else:
                edge_set = (v, u, t)

            if edge_set not in self.sigma_t:
                result = self.compute_sigma_at_one_time(u, v, t)
                # print u,v,t,result
                self.sigma_t[edge_set] = result
                if result > self.eps:
                    tau += 1
            else:
                if self.sigma_t[edge_set] > self.eps:
                    tau += 1

            if tau >= self.tau:
                return tau
        return 0

    def compute_sigma2(self, u, v):
        tau = 0
        if len(self.G.edges[u, v]['t']) < self.tau:
            return 0

        for t in self.G.edges[u, v]['t']:
            if t in self.subgraph:
                pass
            else:
                self.subgraph[t] = nx.subgraph.subgraph(self.path, t, self.theta)
                print(t)
            if self.subgraph[t].compute_sigma(u, v, self.eps) > self.eps:
                tau = tau + 1
            if tau >= self.tau:
                return tau
        return 0

    def SCANA(self, miu, tau, eps):
        self.eps = eps
        self.miu = miu
        self.tau = tau
        self.union_set = []
        #for u in self.G.nodes():
            #### this part has done already in tGraph
            #wt = 0
            #self.G.nodes[u]['l'] = 0
            #wt = 0
            #for nbr, eattr in self.G.adj[u].items():
            #    wt = wt + round(mean(eattr['w']), 2)
            #self.G.nodes[u]['u'] = len(self.G.adj[u])
        ####
        counter = 0
        starttime = datetime.datetime.now()
        """cpu_to_relax = 3
        #pool = ThreadPool(processes=multiprocessing.cpu_count() - cpu_to_relax)
        pool = ThreadPool(processes=1)
        results = pool.starmap(self.SCANA_multiprocess, zip(itertools.repeat(self), self.rank))
        pool.close()
        pool.join()"""
        for u in self.rank:
            counter += 1
            print(counter)
            value = self.check_SCANA_core(u)
            if self.G.nodes[u]['l'] >= self.miu:
                self.add_node_set(u)
                self.cluster_SCANA_core(u)
        endtime = datetime.datetime.now()
        interval = (endtime - starttime).total_seconds()
        print("Runing time of SCANA:" + str(interval))
        self.write_runtime(interval, sys._getframe().f_code.co_name)
        self.sigma = {}
        self.sigma_t = {}
        self.visited_node = []
        file_name = self.path + '.output-' + str(self.eps) + '-' + str(self.tau) + '-' + str(self.miu) + '_SCANA'
        print("Cores output at: " + file_name)
        print(self.union_set)
        file_object = open(file_name, 'w')
        for unit in self.union_set:
            file_object.write(json.dumps(unit))
            file_object.write("\n")
        file_object.close()
        #self.union_set = []

    def cluster_by_cores(self, file_c, flag):
        self.union_set = []
        nodes_set = []

        if flag == 1:
            for line in open(file_c):
                number = int(line.strip('\n'))
                nodes_set.append(number)
        if flag == 2:
            for line in open(file_c):
                new_dict = json.loads(line)
                for item in new_dict:
                    nodes_set.append(item)

        for u in nodes_set:
            self.add_node_set(u)
            for v in set(nodes_set) & set(self.G.adj[u]):
                if v < u:
                    self.add_node_set(v)
                    result = self.compute_sigma(u, v)
                    if result >= self.tau:
                        self.union(u, v)

        cluster_ans = []
        for i in range(len(self.union_set)):
            cluster1 = self.union_set[i][:]
            for u in self.union_set[i]:
                for v in set(self.G.adj[u]):
                    if v <= u:
                        edge_set = (u, v)
                    else:
                        edge_set = (v, u)
                    if edge_set not in self.sigma:
                        result = self.compute_sigma(u, v)
                        self.sigma[edge_set] = result
                    if self.sigma[edge_set] >= self.tau:
                        cluster1.append(v)
            cluster_ans.append(set(cluster1))
        # pprint.pprint(self.union_set)
        file_object = open(file_c + '_cluster', 'w')
        for unit in cluster_ans:
            file_object.write(json.dumps(list(unit)))
            file_object.write("\n")
        file_object.close()

        self.union_set = []
        self.sigma = {}
        self.sigma_t = {}
        self.visited_node = []

    def run(self, filename):
        self.SCANA(self.miu, self.tau, self.eps)

    def cluster(self, filename):
        newname = filename + '.output-' + str(self.eps) + '-' + str(self.tau) + '-' + str(self.miu) + '_SCANA'
        self.cluster_by_cores(newname, 2)

    def write_runtime(self, t, module_name):
        file_object = open('running time', 'a')
        time = {"name": self.path, "eps": self.eps, "tau": self.tau, "miu": self.miu, "time": t,
                "method": module_name}
        file_object.write(json.dumps(time))
        file_object.write("\n")
        file_object.close()


if __name__ == '__main__':
    # get the dataset
    filename = "temporal_df_2020"
    # convert the temporal dataset separated by \t to a networkx. calculate the degree of nodes and find each
    # node's adjacency
    # !call the instructor of the class to initiate members
    # !in this method, these members will be initiated: rank; adj
    # !rank: list of nodes sorted based on degree
    # !adj: dict {node: [list of neighbors sorted based on degree].. }
    G = tGraph(filename)
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    # We can set the number of bins with the `bins` kwarg
    degree_list = []
    for node, deg in G.adj.items():
        degree_list.append(len(deg))
    axs.hist(degree_list, bins=10)
    plt.show()

    # set parameters
    G.eps = 0.5
    G.tau = 1
    G.miu = 500

    print(filename, G.eps, G.tau, G.miu)

    #G.run(filename)
    G.cluster(filename)

    #G.nodes_distribution_by_year()
    #G.degree_distribution_of_nodes_detemporal()
    #G.degree_distribution_of_nodes()

    #G.evauluation(filename, "0.5-3-3500")
    # G.evaluaition_by_year(filename, "0.5-3-5")

    #G.analyse()
