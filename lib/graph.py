import networkx as nx
import json
import copy
import math
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

class Interactome(object):

    def __init__(self, queryTopo):
        self.queryTopo = queryTopo

    def drawGraph(self):
        G=nx.Graph()
        for interaction in self.queryTopo.getEdges(blacklist=None):
            G.add_edge(interaction['lowQuery'], interaction['highQuery'], 
                       lowQueryParam = [lowQueryEval for lowQueryEval in interaction['loQueryEval']] ,
                       highQueryParam = [highQueryEval for highQueryEval in interaction['hiQueryEval']])
        
        # Remove Node with no interactions -- /!\ A REVOIR
        for node in G.node.keys():
            if G.neighbors(node) <= 0:
                G.remove_node(node)
        return G

    def createNeiGraph(self, query, graphParent):
        neighboorGraph = NeighboorGraph(query, graphParent) # /!\ A FUSIONNER AVEC DRAWGRAPH SI ON PEUT plt.
        neighboorGraph.generateNeighboorGraph()
        return neighboorGraph

    def drawNeiGraph(self, neighborsGraph):
        return neighborsGraph.graph

    def filterGraph(self, neighborsGraph, **kwargs):
    
        graph = neighborsGraph.graph

        coverage = 0
        identity = 0
        
        for param in kwargs:
            if param == 'coverage':
                coverage = kwargs['coverage']
            if param == 'identity':
                identity = kwargs['identity']
        
        G = copy.deepcopy(graph)
        
        for edge in G.edge:
            for node in G[edge].keys():
                for i, lowQueryEval in enumerate(G[edge][node]['lowQueryParam']):
                    coverPerCent =  self.coverageCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[5])
                    if float(lowQueryEval[4]) > float(kwargs['evalue']) or float(coverPerCent) < float(coverage):
                        G.adj[edge][node]['highQueryParam'].pop(i)
                        G.adj[edge][node]['lowQueryParam'].pop(i)
                        break
                
                #if node in G[edge] and len(G[edge][node]['highQueryParam']) > 0:
                for i, highQueryEval in enumerate(G[edge][node]['highQueryParam']):
                    coverPerCent =  self.coverageCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[5])

                    if  float(highQueryEval[4]) > float(kwargs['evalue']) or float(coverPerCent) < float(coverage):
                        G.adj[edge][node]['highQueryParam'].pop(i)
                        G.adj[edge][node]['lowQueryParam'].pop(i)
                        break

                # If there is two highQ or lowQ value for the same node
                if len(G.adj[edge][node]['highQueryParam']) == 0 or len(G.adj[edge][node]['lowQueryParam']) == 0:
                    del G.adj[edge][node]
                    
        # Remove Node with no interactions
        for node in G.node.keys():
            if not G.neighbors(node):
                G.remove_node(node)

        neighborsGraph.graph = G
        return neighborsGraph

    def coverageCalculation(self, minSeq, maxSeq, totalSeq):
        coverPerCent =  int(((float(maxSeq) - float(minSeq)) / float(totalSeq)) * 100)
        return coverPerCent


    def drawCurveParam(self, neighborsGraph):
        '''
        for k, v in neighborsGraph.graph.edge.iteritems():
            for value in v:
                print k, value, neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][4],neighborsGraph.graph.edge[k][value]['highQueryParam'][0][4], self.coverageCalculation(neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][0],
                 neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][1], 
                 neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][5])
        '''

        graph = neighborsGraph.graph
        queryCenter = neighborsGraph.queryCenter

        all_evalue = []
        all_coverage = []
        all_identity = []
        
        for query in graph.edge:
            if query.query == queryCenter:
                for neighbor, param in graph.edge[query].iteritems():
                    for low in param['lowQueryParam']:
                        coverPerCent = self.coverageCalculation(low[0], low[1], low[5])
                        all_evalue.append(math.log10(float(low[4])))
                        all_coverage.append(coverPerCent)
                    
                    
                    for high in param['highQueryParam']:
                        coverPerCent = self.coverageCalculation(high[0], high[1], high[5])
                        # Transform to log scale for th evalue
                        all_evalue.append(math.log10(float(high[4])))
                        all_coverage.append(coverPerCent)
                    
        
        sorted_all_evalue_per_cent = sorted([ev for ev in all_evalue])
        sorted_all_coverage = sorted(all_coverage)

        all_evalue = sorted([evalue[1] for evalue in self.densityCumulate(sorted_all_evalue_per_cent)])
        all_coverage = sorted([evalue[1] for evalue in self.densityCumulate(sorted_all_coverage)])

        cumulate_node_evalue = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_evalue_per_cent)])
        cumulate_node_coverage = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_coverage)])

        #Legend
        red_patch = mpatches.Patch(color='red', label='eValue')
        green_patch = mpatches.Patch(color='green', label='Coverage')
        
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111)
        ax1.legend(handles=[red_patch])
        ax1.plot(all_evalue, cumulate_node_evalue, color = 'r', linestyle='--', marker='o')
        plt.show()

        fig2 = plt.figure(2)
        ax2 = fig2.add_subplot(111)
        ax2.legend(handles=[green_patch])
        ax2.plot(all_coverage, cumulate_node_coverage, linestyle='--', marker='o', color='g')
        plt.show()

    def densityCumulate(self, listParam):
        dataSet = set()
        for param in listParam:
            dataSet.add((len([nodeFiltered for nodeFiltered in listParam if nodeFiltered <= param]), param))
        return dataSet

    def drawNeiTopo(self, neighbors_dict):
        print "Liste des 1ers voisine:\n"
        for node in neighbors_dict:
            print ', '.join([neighbor.query for neighbor in neighbors_dict[node]])

    def serializeGraph(self, queryCenter, graphEdge, path):
    
        jsonStruct = {"Queries" : {}}
        for query, nodes in graphEdge.iteritems():
            if query.query == queryCenter:
                jsonStruct["Queries"][query.query] = {}
                for node, param in nodes.iteritems():
                    jsonStruct["Queries"][query.query][node.query] = {"lowQueryParam" : [low for low in param["lowQueryParam"]],
                                                                      "highQueryParam" : [high for high in param["highQueryParam"]]}
                
        json.dump(jsonStruct, file(path, 'w'))

    def deserializeGraph(self, beanPath):
        G=nx.Graph()
        with open (beanPath, 'r') as file:
            data = json.load(file)
            for query in data['Queries']:
                for neighbor, param in data['Queries'][query].iteritems():
                        # Creation de edges
                    G.add_edge(self.queryTopo.dictQuery[query][0], self.queryTopo.dictQuery[neighbor][0], 
                               lowQueryParam = [low for low in param['lowQueryParam']],
                               highQueryParam = [high for high in param['highQueryParam']])
        return G

class NeighboorGraph(object):
    def __init__(self, queryCenter, graphParent):
        self.graphParent = graphParent
        self.queryCenter = queryCenter
        self.graph = None

    def generateNeighboorGraph(self):
        G=nx.Graph()
        
        # Declaration des variables
        neighborParam = []
        queryParam = []
        
        for query in self.graphParent.edge:
            if query.query == self.queryCenter:
                for neighbor, param in self.graphParent.edge[query].iteritems():
                    # Creation de edges
                    G.add_edge(query, neighbor, 
                        lowQueryParam = param['lowQueryParam'],
                        highQueryParam = param['highQueryParam'])
        self.graph = G
        return self.graph