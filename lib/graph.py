import networkx as nx
import json
import copy
import math
import numpy as np
import random
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

    def filterGraph(self, neighborsGraph, **kwargs):
    
        graph = neighborsGraph.graphData

        coverage = 1
        identity = 1
        similarity = 1
        
        for param in kwargs:
            if param == 'coverage':
                coverage = kwargs['coverage']
            if param == 'identity':
                identity = kwargs['identity']
            if param == 'similarity':
                identity = kwargs['similarity']
        
        G = copy.deepcopy(graph)
        
        for edge in G.edge:
            for node in G[edge].keys():
                for i, lowQueryEval in enumerate(G[edge][node]['lowQueryParam']):
                    
                    coverPerCent =  self.coverageCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[5])
                    idPerCent = self.idCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[3])
                    simiPerCent = self.simiCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[2])
                    
                    if float(lowQueryEval[4]) > float(kwargs['evalue']) or float(coverPerCent) < float(coverage) or idPerCent < identity or simiPerCent < similarity:
                        G.adj[edge][node]['highQueryParam'].pop(i)
                        G.adj[edge][node]['lowQueryParam'].pop(i)
                        break
                
                #if node in G[edge] and len(G[edge][node]['highQueryParam']) > 0:
                for i, highQueryEval in enumerate(G[edge][node]['highQueryParam']):
                    
                    coverPerCent =  self.coverageCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[5])
                    idPerCent = self.idCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[3])
                    simiPerCent = self.simiCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[2])
                    
                    if  float(highQueryEval[4]) > float(kwargs['evalue']) or float(coverPerCent) < float(coverage) or idPerCent < identity or simiPerCent < similarity:
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

        neighborsGraph.graphData = G
        return neighborsGraph

    def coverageCalculation(self, minSeq, maxSeq, totalSeq):
        coverPerCent =  int(((float(maxSeq) - float(minSeq)) / float(totalSeq)) * 100)
        return coverPerCent

    def idCalculation(self, minSeq, maxSeq, idValue):
        idPerCent =  int((float(idValue) / (float(maxSeq) - float(minSeq))) * 100)
        return idPerCent

    def simiCalculation(self, minSeq, maxSeq, simiValue):
        simiCalculation =  int((float(simiValue) / (float(maxSeq) - float(minSeq))) * 100)
        return simiCalculation

    def drawCurveParam(self, neighborsGraph):
        '''
        for k, v in neighborsGraph.graph.edge.iteritems():
            for value in v:
                print k, value, neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][4],neighborsGraph.graph.edge[k][value]['highQueryParam'][0][4], self.coverageCalculation(neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][0],
                 neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][1], 
                 neighborsGraph.graph.edge[k][value]['lowQueryParam'][0][5])
        '''

        graph = neighborsGraph.graphData
        queryCenter = neighborsGraph.queryCenter

        all_evalue = []
        all_coverage = []
        all_identity = []
        all_similarity = []
        
        for query in graph.edge:
            if query.query == queryCenter:
                for neighbor, param in graph.edge[query].iteritems():
                    for low in param['lowQueryParam']:
                        
                        coverPerCent = self.coverageCalculation(low[0], low[1], low[5])
                        idPerCent = self.idCalculation(low[0], low[1], low[3])
                        simiPerCent = self.simiCalculation(low[0], low[1], low[2])

                        all_evalue.append(math.log10(float(low[4])))
                        all_coverage.append(coverPerCent)
                        all_identity.append(idPerCent)
                        all_similarity.append(simiPerCent)
                    
                    
                    for high in param['highQueryParam']:
                        coverPerCent = self.coverageCalculation(high[0], high[1], high[5])
                        idPerCent = self.idCalculation(high[0], high[1], high[3])
                        simiPerCent = self.simiCalculation(high[0], high[1], high[2])
                        
                        # Transform to log scale for th evalue
                        all_evalue.append(math.log10(float(high[4])))
                        all_coverage.append(coverPerCent)
                        all_identity.append(idPerCent)
                        all_similarity.append(simiPerCent)
                    
        
        sorted_all_evalue_per_cent = sorted([ev for ev in all_evalue])
        sorted_all_coverage = sorted(all_coverage)
        sorted_all_identity = sorted(all_identity)
        sorted_all_similarity = sorted(all_similarity)

        all_evalue = sorted([evalue[1] for evalue in self.densityCumulate(sorted_all_evalue_per_cent)])
        all_coverage = sorted([evalue[1] for evalue in self.densityCumulate(sorted_all_coverage)])
        all_identity = sorted([evalue[1] for evalue in self.densityCumulate(sorted_all_identity)])
        all_similarity = sorted([evalue[1] for evalue in self.densityCumulate(sorted_all_similarity)])

        cumulate_node_evalue = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_evalue_per_cent)])
        cumulate_node_coverage = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_coverage)])
        cumulate_node_identity = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_identity)])
        cumulate_node_similarity = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_similarity)])

        #Legend
        red_patch = mpatches.Patch(color='red', label='eValue')
        green_patch = mpatches.Patch(color='green', label='Coverage')
        blue_patch = mpatches.Patch(color='blue', label='Identity')
        yellow_patch = mpatches.Patch(color='yellow', label='Similarity')
        
        fig1 = plt.figure(random.randint(0, 1000))
        ax1 = fig1.add_subplot(111)
        ax1.legend(handles=[red_patch])
        ax1.plot(all_evalue, cumulate_node_evalue, color = 'r', linestyle='--', marker='o')
        plt.show()

        fig2 = plt.figure(random.randint(0, 1000))
        ax2 = fig2.add_subplot(111)
        ax2.legend(handles=[green_patch, blue_patch, yellow_patch])
        ax2.plot(all_coverage, cumulate_node_coverage, linestyle='--', marker='o', color='g')
        ax2.plot(all_identity, cumulate_node_identity, linestyle='--', marker='o', color='b')
        ax2.plot(all_similarity, cumulate_node_similarity, linestyle='--', marker='o', color='y')
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

class NeighboorGraph(object):
    def __init__(self, queryCenter, graphParent):
        self.graphParent = graphParent
        self.queryCenter = queryCenter
        self.neighborsList = []
        self.graphData = None

    def generateNeighboorGraph(self, **kwargs):
        G=nx.Graph()
        fig1 = plt.figure(random.randint(0, 1000))
        ax1 = fig1.add_subplot(111)
        
        # Declaration des variables
        neighborParam = []
        queryParam = []
        for query in self.graphParent.edge:
            if query.query == self.queryCenter:
                for neighbor, param in self.graphParent.edge[query].iteritems():
                    
                    # Append all neighboors ID
                    self.neighborsList.append(neighbor.query)

                    # Creation de edges
                    G.add_edge(query, neighbor, 
                        lowQueryParam = param['lowQueryParam'],
                        highQueryParam = param['highQueryParam'])

        # Build 1st Neighboor edges
        if "isConnected" in kwargs:
            if kwargs['isConnected']:
                for edge in self.retrieveNeighboorsEdges(self.graphParent.edge, self.neighborsList):
                    G.add_edge(edge[0], edge[1], 
                                lowQueryParam = edge[2],
                                highQueryParam = edge[3])

        self.graphData = G
        return

    def retrieveNeighboorsEdges(self, fullGraphEdge, listNeighboors):
        
        listNeighboorsEdges = []

        for query in self.graphParent.edge:
            if query.query in listNeighboors:
                for neighbor, param in self.graphParent.edge[query].iteritems():
                    if neighbor.query in listNeighboors:
                        yield(query, neighbor, param['lowQueryParam'], param['highQueryParam'])

class dataOperations(object):
    def __init__(self, queryData):
        # queryData link the the query string to a query memory adress
        self.queryData = queryData

    def serializeGraph(self, neighboorGraph, path):

        graphEdge = neighboorGraph.graphData.edge
        queryList = neighboorGraph.neighborsList
        queryCenter = neighboorGraph.queryCenter

        jsonStruct = {"Queries" : {}}

        for query, nodes in graphEdge.iteritems():
            if query.query in queryList or query.query in queryCenter:
                jsonStruct["Queries"][query.query] = {}
                for node, param in nodes.iteritems():
                    jsonStruct["Queries"][query.query][node.query] = {"lowQueryParam" : [low for low in param["lowQueryParam"]],
                                                                      "highQueryParam" : [high for high in param["highQueryParam"]]}
        json.dump(jsonStruct, file(path, 'w'))

    
    def deserializeGraph(self, beanPath):
        G=nx.Graph()
        fig1 = plt.figure(random.randint(0, 1000))
        ax1 = fig1.add_subplot(111)

        with open (beanPath, 'r') as file:
            data = json.load(file)
            for query in data['Queries']:
                for neighbor, param in data['Queries'][query].iteritems():
                        # Creation de edges
                    G.add_edge(self.queryData.dictQuery[query][0], self.queryData.dictQuery[neighbor][0], 
                               lowQueryParam = [low for low in param['lowQueryParam']],
                               highQueryParam = [high for high in param['highQueryParam']])

        nx.draw_networkx(G, with_labels = True)
        plt.show

        return 















