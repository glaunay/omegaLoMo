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
    
        coverage = 0
        identity = 0
        
        for param in kwargs:
            if param == 'coverage':
                coverage = kwargs['coverage']
            if param == 'identity':
                identity = kwargs['identity']
        
        G = copy.deepcopy(neighborsGraph)
        
        for edge in G.edge:
            for node in G[edge].keys():
                for i, lowQueryEval in enumerate(G[edge][node]['lowQueryParam']):
                        
                    coverPerCent = (int(lowQueryEval[3]) * 100) / (int(lowQueryEval[1]) - int(lowQueryEval[0])) 
                   
                    if float(kwargs['evalue']) > float(lowQueryEval[4]) and float(coverage) > float(coverPerCent): 
                        G.adj[edge][node]['highQueryParam'].pop(i)
                        G.adj[edge][node]['lowQueryParam'].pop(i)

                        break
                
                if node in G[edge] and len(G[edge][node]['highQueryParam']) > 0:
                    for i, highQueryEval in enumerate(G[edge][node]['highQueryParam']):
                        
                        coverPerCent = (int(highQueryEval[3]) * 100) / (int(highQueryEval[1]) - int(highQueryEval[0]))

                        if float(kwargs['evalue']) > float(highQueryEval[4]) and float(coverage) > float(coverPerCent):
                            G.adj[edge][node]['highQueryParam'].pop(i)
                            G.adj[edge][node]['lowQueryParam'].pop(i)
                            break

                # If there is two highQ or lowQ value for the same node
                if len(G.adj[edge][node]['highQueryParam']) <= 0 or len(G.adj[edge][node]['lowQueryParam']) <= 0:
                    del G.adj[edge][node]
                    
        # Remove Node with no interactions
        for node in G.node.keys():
            if not G.neighbors(node):
                G.remove_node(node)

        return G

    def drawCurveParam(self, neighborsGraph):
        
        graph = neighborsGraph.graph
        queryCenter = neighborsGraph.queryCenter

        all_evalue = set()
        all_coverage = set()
        all_identity = []
        
        for query in graph.edge:
            if query.query == queryCenter:
                for neighbor, param in graph.edge[query].iteritems():
                    for low in param['lowQueryParam']:
                       
                        coverPerCent = (int(low[3]) * 100) / (int(low[1]) - int(low[0])) 
                        all_evalue.add(math.log10(float(low[4])))
                        all_coverage.add(coverPerCent)
                    
                    for high in param['highQueryParam']:
                        
                        coverPerCent = (int(high[3]) * 100) / (int(high[1]) - int(high[0])) 
                        # Transform to log scale for th evalue
                        all_evalue.add(math.log10(float(high[4])))
                        all_coverage.add(coverPerCent)
                        
        f1 = plt.figure()
        ax1 = f1.add_subplot(111)
        
        #Legend
        red_patch = mpatches.Patch(color='red', label='eValue')
        green_patch = mpatches.Patch(color='green', label='Coverage')
        
        ax1.legend(handles=[red_patch, green_patch])
        
        xEval = np.arange(0, len(all_evalue), 1)
        xCov = np.arange(0, len(all_coverage), 1)
        
        all_evalue = sorted(all_evalue)
        all_coverage = sorted(all_coverage)
        
        ax1.plot(xEval, all_evalue, color = 'r', linestyle='--', marker='o')
        # Il se peut que le nombre de coverage sois inferieur au nombre de evalue car ecrasement de valeurs
        ax1.plot(xCov, all_coverage, linestyle='--', marker='o', color='g')
        return ax1

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
        return G