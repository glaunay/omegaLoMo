import networkx as nx
import json
import copy
import math
import numpy as np
import random
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt



'''
A utility class to facilitate the visualisation of graph statistic in tabular format under jupyter
'''

class NodeView(object):
    def __init__(self, G):
        self.G = G
    def _repr_html_(self):
        htmlString = '<table><thead><th>Node</th><th>Degree</th></thead><tbody>'
        d = nx.degree(self.G)
        for n in sorted( [ { 'nodeRef' : k, 'degree' : d[k] } for k in d ], key=lambda k: k['degree'], reverse=True):
            color = 'brown' if self.G.node[n['nodeRef']]['group'] == 1 else 'black'

            htmlString += '<tr><td><a href="http://www.uniprot.org/uniprot/' + str(n['nodeRef']) + '" target="blank" style="color:' + color + '">' + str(n['nodeRef']) + '</td>'
            htmlString += '<td>' + str(n['degree']) + '</td></tr>'

        htmlString += '</tbody></table>'

        return htmlString


class Interactome(object):

    def __init__(self, queryTopo):
        self.queryTopo = queryTopo
        self.singleNode = {}
        self._Gtotal = None
        self.verbose=False
        self._nodeView = None

    def coverageCalculation(self, minSeq, maxSeq, totalSeq):
        coverPerCent =  int((((float(maxSeq) - float(minSeq)) + 1) / float(totalSeq)) * 100)
        return coverPerCent

    def idCalculation(self, minSeq, maxSeq, idValue):
        idPerCent =  int((float(idValue) / (float(maxSeq) - float(minSeq))) * 100)
        return idPerCent

    def simiCalculation(self, minSeq, maxSeq, simiValue):
        simiCalculation =  int((float(simiValue) / (float(maxSeq) - float(minSeq))) * 100)
        return simiCalculation

    @property
    def Gtotal(self):
        if not self._Gtotal:
            self._Gtotal = self.createGraph()

        return self._Gtotal

    # Returns a networkX graph object, applying optional filtering to the global graph
    # We have to check for moification to the global graph through references
    def graph(self, **kwargs):
        G = self.Gtotal.copy()

        if 'verbose' in kwargs:
            self.verbose = kwargs['verbose']

       # return G

        seedNodes = kwargs['seedNodes'] if 'seedNodes' in kwargs else None
        radius = kwargs['radius'] if 'radius' in kwargs else 1
        coverage =  kwargs['coverage'] if 'coverage' in kwargs else 1
        simPct =  kwargs['simPct'] if 'simPct' in kwargs else None

# Extracting subnetwork around seeds and merge them
        if seedNodes:
            Gall = []
            #totalNodes = set()
            #scoreBoard = {}
            for n in G.nodes():
                if hash(n) in [ hash(ns) for ns in seedNodes ]:
                    G.node[n]['group'] = 1
                    Gtmp = nx.ego_graph(G, n, radius)
                    Gall.append(Gtmp)
                    if self.verbose:
                        print str(len(Gtmp.nodes())) + " nodes around " + str(n)
            G = Gall[0]
            for i in range(1,len(Gall)):
                G = nx.compose(G, Gall[i])

# Graph utility functions
        def trimEdgeData(G, edge, coverageTreshold=0.00, simTreshold=0.00):
            e=G.get_edge_data(*edge)
            remove_indices = []
            for i, (leftNodeParam, rightNodeParam) in enumerate(zip(e['lowQueryParam'], e['highQueryParam'])):
                leftCoverPerCent =  self.coverageCalculation(leftNodeParam[0], leftNodeParam[1], leftNodeParam[8])
                rightCoverPerCent =  self.coverageCalculation(rightNodeParam[0], rightNodeParam[1], rightNodeParam[8])
                leftSimPerCent =  self.simiCalculation(leftNodeParam[0], leftNodeParam[1], leftNodeParam[2])
                rightSimPerCent =  self.simiCalculation(rightNodeParam[0], rightNodeParam[1], rightNodeParam[2])
                if min(leftCoverPerCent, rightCoverPerCent) < coverageTreshold or min(leftSimPerCent, rightSimPerCent) < simTreshold:
                    remove_indices.append(i)


            if not remove_indices:
                return

            if self.verbose:
                print "removing following edge data indices : " + str(remove_indices)

            e['lowQueryParam'] = [ d for i,d in enumerate (e['lowQueryParam']) if i not in remove_indices ]
            e['highQueryParam'] = [ d for i,d in enumerate (e['highQueryParam']) if i not in remove_indices ]

            return

        def edgeDepth(G, e):
            d=G.get_edge_data(*e)
            return len(d['lowQueryParam'])

        if 'coverage' in kwargs or 'simPct' in kwargs:
            totalEdges = 0
            nodePairToPrune = []

            for e in G.edges():
                totalEdges += 1
                if self.verbose:
                    print 'intial depth of ' + str(e) + ' is ' + str(edgeDepth(G, e))
                trimEdgeData(G, e, coverage, simPct)
                if self.verbose:
                    print 'trimmed depth of ' + str(e) + ' is ' + str(edgeDepth(G, e))
                if edgeDepth(G, e) == 0:
                    if self.verbose:
                        print (str(e) + ' is an empty edge, registered for pruning')
                    nodePairToPrune.append(e[:2])

                    G.remove_edge(*e[:2]) # unpacks e from an edge tuple

        print "Total number of pruned edges is " + str(len(nodePairToPrune)) + ' / ' + str(totalEdges)

        #remove isolates
        G.remove_nodes_from(nx.isolates(G))
        #remove self connected only
        toDel = []
        for n in G.nodes():
            if G.degree(n) > 1:
                continue;
            if G.neighbors(n)[0] == n:
                print n + " is self connected only"
                toDel.append(n)
        #remove self-connected only
        G.remove_nodes_from(toDel)

        return G


    def createGraph(self):
        G=nx.Graph()
        for interaction in self.queryTopo.getEdges(blacklist=None):

            G.add_edge(interaction['lowQuery'], interaction['highQuery'],

                       lowQueryParam = [lowQueryEval for lowQueryEval in interaction['loQueryEval']] ,
                       highQueryParam = [highQueryEval for highQueryEval in interaction['hiQueryEval']])

        # Remove Node with no interactions -- /!\ A REVOIR
        for node in G.node.keys():
            if len(G.neighbors(node)) <= 1:
                self.singleNode[node.query] = node
                G.remove_node(node)

        nx.set_node_attributes(G, 'group', 0)
        #print G.nodes(data=True)

        return G



'''

    def createNeiGraph(self, query, graphParent, draw):
        neighboorGraph = NeighboorGraph(query, graphParent) # /!\ A FUSIONNER AVEC DRAWGRAPH SI ON PEUT plt.
        neighboorGraph.generateNeighboorGraph(draw)
        return neighboorGraph

    def filterGraph(self, neighborsGraph, **kwargs):

        G=nx.Graph()
        fig1 = plt.figure(random.randint(0, 1000))
        ax1 = fig1.add_subplot(111)

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
                similarity = kwargs['similarity']

        G = copy.deepcopy(graph)

        for edge in G.edge:
            for node in G[edge].keys():
                for i, lowQueryEval in enumerate(G[edge][node]['lowQueryParam']):

                    coverPerCent =  self.coverageCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[8])
                    idPerCent = self.idCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[3])
                    simiPerCent = self.simiCalculation(lowQueryEval[0], lowQueryEval[1], lowQueryEval[2])

#                    print str(coverPerCent) + " " + str(idPerCent) + " " + str(simiPerCent)

                    if float(lowQueryEval[4]) > float(kwargs['evalue']) or float(coverPerCent) < float(coverage) or float(idPerCent) < float(identity) or float(simiPerCent) < float(similarity):
                        G.adj[edge][node]['highQueryParam'].pop(i)
                        G.adj[edge][node]['lowQueryParam'].pop(i)
                        break

                #if node in G[edge] and len(G[edge][node]['highQueryParam']) > 0:
                for i, highQueryEval in enumerate(G[edge][node]['highQueryParam']):

                    coverPerCent =  self.coverageCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[8])
                    idPerCent = self.idCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[3])
                    simiPerCent = self.simiCalculation(highQueryEval[0], highQueryEval[1], highQueryEval[2])

                    if  float(highQueryEval[4]) > float(kwargs['evalue']) or float(coverPerCent) < float(coverage) or float(idPerCent) < float(identity) or float(simiPerCent) < float(similarity):
                        G.adj[edge][node]['highQueryParam'].pop(i)
                        G.adj[edge][node]['lowQueryParam'].pop(i)
                        break

        for edge in G.adj.keys():
            if len(G.adj[edge]) > 0:
                for node in G.adj[edge].keys():
                    if len(G.adj[edge][node]['highQueryParam']) == 0 or len(G.adj[edge][node]['lowQueryParam']) == 0 :
                        del G.adj[edge][node]

            else:
                del G.adj[edge]

        neiglist = []
        for node in G.node.keys():
            if not G.neighbors(node):
                G.remove_node(node)
            else:
                neiglist.append(node.query)


        #print G.edge, G.node
        print ', '.join([fn for fn in neiglist])
        return G#, neighborsGraph
'''

'''
    def drawCurveParam(self, neighborsGraph):

        graph = neighborsGraph.graphData
        queryCenter = neighborsGraph.queryCenter

        all_evalue = []
        all_coverage = []
        all_identity = []
        all_similarity = []
        edge_minValue = []

        plotlow = []
        plothig = []

        for query in graph.edge:
            if query.query == queryCenter:
                for neighbor, param in graph.edge[query].iteritems():

                    #print neighbor, param, query
                    for c, i in enumerate(param['lowQueryParam']):
                        if float(param['lowQueryParam'][c][4]) <= 0 or float(param['highQueryParam'][c][4]) <= 0 :
                            #print query, neighbor, param
                            pass
                        else :
                            plotlow.append(math.log10(float(param['lowQueryParam'][c][4])))
                            plothig.append(math.log10(float(param['highQueryParam'][c][4])))
                            if param['lowQueryParam'][c][4] < param['highQueryParam'][c][4] :
                                #print param['lowQueryParam'][c][4]
                                edge_minValue.append(math.log10(float(param['lowQueryParam'][c][4])))
                            else :
                                #print param['highQueryParam'][c][4]
                                edge_minValue.append(math.log10(float(param['highQueryParam'][c][4])))


                            for low in param['lowQueryParam']:

                                coverPerCent = self.coverageCalculation(low[0], low[1], low[8])
                                #if coverPerCent > 100:
                                    #print query.query, query.template, neighbor.query, neighbor.template, coverPerCent, param
                                idPerCent = self.idCalculation(low[0], low[1], low[3])
                                simiPerCent = self.simiCalculation(low[0], low[1], low[2])
                                print simiPerCent

                                all_evalue.append(math.log10(float(low[4])))
                                all_coverage.append(coverPerCent)
                                all_identity.append(idPerCent)
                                all_similarity.append(simiPerCent)


                            for high in param['highQueryParam']:
                                coverPerCent = self.coverageCalculation(high[0], high[1], high[8])
                                #if coverPerCent > 100:
                                #    print query.query, query.template, neighbor.query, neighbor.template, coverPerCent
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
        cumulate_node_coverage = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_coverage)], reverse=True)
        cumulate_node_identity = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_identity)], reverse=True)
        cumulate_node_similarity = sorted([evalue[0] for evalue in self.densityCumulate(sorted_all_similarity)], reverse=True)

        #Legend
        red_patch = mpatches.Patch(color='red', label='eValue')
        green_patch = mpatches.Patch(color='green', label='Coverage')
        blue_patch = mpatches.Patch(color='blue', label='Identity')
        yellow_patch = mpatches.Patch(color='yellow', label='Similarity')

        if not len(graph.edge) == 0 :

            fig1 = plt.figure(random.randint(0, 1000))
            ax1 = fig1.add_subplot(111)
            ax1.legend(handles=[red_patch], bbox_to_anchor=(1, 1), loc=2)
            ax1.plot(all_evalue, cumulate_node_evalue, color = 'r', linestyle='--', marker='o')
            plt.show()

            #hist eV
            plt.hist([all_evalue], bins=20 , color=['red'], label=['eV'])
            plt.legend()
            plt.show()
            #hist eV edges
            plt.hist([edge_minValue], bins=20 , color=['red'], label=['eVEdge'])
            plt.legend()
            plt.show()

            fig2 = plt.figure(random.randint(0, 1000))
            ax2 = fig2.add_subplot(111)
            ax2.legend(handles=[green_patch, blue_patch, yellow_patch],bbox_to_anchor=(1, 1), loc=2)

            ax2.plot(all_coverage, cumulate_node_coverage, linestyle='--', marker='o', color='g')
            ax2.plot(all_identity, cumulate_node_identity, linestyle='--', marker='o', color='b')
            ax2.plot(all_similarity, cumulate_node_similarity, linestyle='--', marker='o', color='y')
            ##LOIC
            plt.gca().invert_xaxis()
            ###
            plt.show()


            #hist cover/simi/ident
            plt.hist([all_coverage, all_similarity, all_identity], bins=7, color=['green', 'yellow', 'blue'],
                label=['Coverage', 'Similarity', 'Identity'])
            plt.legend()
            plt.gca().invert_xaxis()
            plt.show()

            plt.plot(plotlow, plothig, 'ro')
            plt.axis([-140, 20 , -140, 20])
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
'''

class NeighboorGraph(object):
    def __init__(self, queryCenter, graphParent):
        self.graphParent = graphParent
        self.queryCenter = queryCenter
        self.neighborsList = []
        self.graphData = None

    def generateNeighboorGraph(self, draw = True):
        G=nx.Graph()
        if draw :
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

        self.graphData = G
        return

    def connectNeighboor(self, neighboorGraph):
        copyNeighboorGraph = copy.deepcopy(neighboorGraph)
        G = copyNeighboorGraph.graphData
        #fig1 = plt.figure(random.randint(0, 1000))
        #ax1 = fig1.add_subplot(111)
        # Build 1st Neighboor edges
        self._generateNeighboorList(copyNeighboorGraph)
        for edge in self.retrieveNeighboorsEdges(copyNeighboorGraph.graphParent.edge, copyNeighboorGraph.neighborsList):
            G.add_edge(edge[0], edge[1],
                lowQueryParam = edge[2],
                highQueryParam = edge[3])

        nx.draw_networkx(copyNeighboorGraph.graphData, with_labels = True)
        plt.show()

        return copyNeighboorGraph


    def retrieveNeighboorsEdges(self, fullGraphEdge, listNeighboors):

        listNeighboorsEdges = []

        for query in self.graphParent.edge:
            if query.query in listNeighboors:
                for neighbor, param in self.graphParent.edge[query].iteritems():
                    if neighbor.query in listNeighboors:
                        yield(query, neighbor, param['lowQueryParam'], param['highQueryParam'])

    def getEdgeInformation(self, idNeighboor):

        dict_edge = self.graphData.edge

        for query in dict_edge:
            if query.query == self.queryCenter:
                for neighboor, param in dict_edge[query].iteritems():
                    if neighboor.query == idNeighboor:
                        #print query, neighboor, param
                        if param['lowQueryParam'][0][4] < param['highQueryParam'][0][4]:
                            return param['lowQueryParam'][0]
                        #    print query, neighboor, param['lowQueryParam'][0][4]
                        else:
                        #    print query, neighboor, param['highQueryParam'][0][4]
                            return param['highQueryParam'][0]

    def _generateNeighboorList(self, neighboorGraph):

        neighboorGraph.neighborsList[:] = []

        for query in neighboorGraph.graphData.edge:
            if query.query == neighboorGraph.queryCenter:
                for neighbor, param in neighboorGraph.graphData.edge[query].iteritems():
                    # Append all neighboors ID
                    neighboorGraph.neighborsList.append(neighbor.query)


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
        plt.show()

        return















