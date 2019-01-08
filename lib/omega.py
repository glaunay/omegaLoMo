import xml.etree.ElementTree as ET
import sys
import json
import numpy as np
from os import listdir
from os.path import isfile, join
import os
import fnmatch
import gzip
from hashlib import md5 as md5
import networkx as nx

import pyproteins.container.Core as co
'''
.....
GL revised version of omega Matrix/Vector implemntations
....
'''

''' OUTPUT 
    (mitabTopologyTree, 
    { link : [ {
        { "source" : NODE_i, "target": NODE_j , "data" : { "lowQueryParam" : [[]], "highQueryParam"}}
    }]
    ,
    nodes : [ {
        "id" : ,
        "group" : 
        "val" 
    }]
    
    )


(Q8DR57, #Node uniprot Identifier on the "lowQuery" side of the edge
 Q9EUQ7, #Node uniprot Identifier on the "highQuery" side of the edge
 {'lowQueryParam': [
     # Each homolog relationship supporting interaction inference on the "low" side of the edge
     [u'1', u'105', u'66', u'49', u'3.02566e-17', u'98', u'3', u'98', u'109', 'P0CI74']
     ],
  'highQueryParam': [
      # Each homolog relationship supporting interaction inference on the "high" side of the edge
      [u'3', u'230', u'138', u'84', u'7.50339e-38', u'251', u'6', u'242', u'242', 'P35154']
      ]
})
'''




class homologTree(object):
    def __init__(self, data=None, fileName=None):
        if fileName:
            with open(fileName) as f:
                self.data = json.load(f)
        else:
            raise ValueError("No input")
    '''
    parentID -> Dict
    Where, Dict = 
    {
        childID : [homologyInfo[] + parentID ],

    }
    Two level dict:
    The top key is the uniclust ID
    The 2nd level key is the R6 identifier
    NB: The deep value is the homology vector w/ 
        the following content
     [ R6_sequence Length,
      R6_hsp start position,
      R6_hsp stop position,
      PSQ_sequence Length,
      PSQ_hsp start position,
      PSQ_hsp stop position,
      HSP positive match number,
      HSP identical match number,
      HSP eValue
      ]
    '''
        
    # More explicitely Given a psq id
    # returns a dict W/ keys as R6  identifiers and value a tuple
    # where the original  homology vectors is prepend with the psq id
    # EX:
    #"A0A142I5B9" : {"P0A4D8" : [ ["350","47","347","399","49","352","5.06522e-30"] ], .. }
    # "P0A4D8" : ("A0A142I5B9","47","347","99","49","5.06522e-30"),
    # "Q8DQH1" :("A0A142I5B9","47","363","108","49","1.21959e-20")
    # }

    def getChildrenData(self, psqId):
        if not psqId in self.data:
            return {}
        d = {}
        for homologKey in self.data[psqId]:
            d[homologKey] = [psqId]
            hVector = self.data[psqId][homologKey][0] ## Take 1st HSP if any more
            for e in hVector:
                d[homologKey].append(e)
        return d


## Needs documentation and notebook API use case

class omegaTopology(object):
    def __init__(self, mitabTopologyObject, homologySupportFile):
        self.hData = homologTree(fileName=homologySupportFile)
        self.baseTopology = mitabTopologyObject
        self.adjTree = co.mdTree(append=False)
    
    def jupyterNodeView(self): 
        return NodeView(self._G)

    def resetNodeVisbility(self):
        pass

    def prune(self, *seeds):
        self._G = self._g()
        nx.set_node_attributes(self._G, 0, 'group')

        if not seeds:
            print("No seed to prune")
            #return self._G
        else:
            seedSet = []
            otherSet = []
            for n in self._G.nodes(data=True):
                if n[0] in seeds:
                    n[1]['group'] = 1
                    seedSet.append(n[0])
                else:
                    otherSet.append(n[0])
            
            for node in otherSet:
                connected = False
                for seed in seedSet:
                    if nx.has_path(self._G, node, seed):
                        connected = True
                        break
                if not connected:
                    #print('Removing ', node)
                    self._G.remove_node(node)
                    self._hideNode(node)

        degrees = self._G.degree()
        for n,d in self._G.nodes(data=True):
            d['val'] = degrees[n]

        return self._G
    
    def dump(self):
        nodes = [ { "id" : str(n), "group" : d["group"], "val" : d["val"]}  for n,d in self._G.nodes(data=True) ]        
        links = [ {"source" : x, "target" : y , "data" : d } for x, y, d in self.iterVisible() ]
        return { "nodes" : nodes , "links": links }

    def _g(self):
        G=nx.Graph()
        for n1, n2, edgeData in self.iterVisible():
            #print (n1,n2)
            G.add_edge( n1, n2, data=edgeData )
        return G
        
    ## List all nodes as primary key and set of templates as values
    @property
    def nodes(self):
        nodes = {}
        #for n1, n2, e in self:
        #    if e.isEmpty:
        #        continue
        for n1, n2, e in self.iterVisible():   
            templates = e.templates
            if n1 not in nodes:
                nodes[n1] = set()
            nodes[n1].update(templates[0])
            
            #if n1 != n2: Not needed
            if n2 not in nodes:
                nodes[n2] = set()
            nodes[n2].update(templates[1])

        return nodes  

    def _hideNode(self, node): # We hide all egdes connecting node
        for partner, edge in self.adjTree[node].items():
            edge.visible = False
        
# Build edge specifications
# Edges carry homology relationship
    def buildEdges(self):
        nbBase = 0
        nbModel = 0
        for baseIdA, baseIdB, datum in self.baseTopology:
            nbBase += 1           
            dataNewA = self.hData.getChildrenData(baseIdA)
            dataNewB = self.hData.getChildrenData(baseIdB)
            self.addEdgeSet(dataNewA, dataNewB)

                
        print(self.edgeNumber, " interactions unpacked from ", nbBase)
       # for x in  self.adjTree:

    # B/C Effectively change the toplogy, 
    # edges are set back to visible=True
    # seed pruning will have to be made again
    def trimEdges(self, simPct = 0.0, idPct = 0.0, cvPct = 0.0):
        nDel = 0
        nTot = 0
        for node1, node2, HoParameterSetObj in self:
            nTot += 1
            HoParameterSetObj.trim(simPct=simPct, idPct=idPct, cvPct=cvPct)
            if HoParameterSetObj.isEmpty:
                nDel += 1
        print(nDel, ' interactions trimmed from total ', nTot)

    def iterVisible(self):
        for x,y,e in self:
            if not e.isEmpty and e.visible:
                yield (x,y,e)

    def __iter__(self):
        for k1,k2, edgeDatum in self.adjTree:
            yield ( k1, k2, edgeDatum)

    @property
    def edgeNumber(self):
        return len([ e for x,y,e in self.iterVisible() ])

    def __repr__(self):
        s = []
        for n1, n2, e in self.iterVisible():
            s.append( (n1, n2, e) ) 
        
        return str(s)

    def getEdgeSet(self, *args):
        if not args:
            return self.iterVisible()
        if len(args) > 2:
            raise ValueError('expected 0, 1 or 2 node got', len(args))


    ## Storing adjacency and mitab evidences in a two level dict
    # using md5 sort sort keys, so that JS implementation can rely on md5 too   
    def addEdgeSet(self, dataNewA, dataNewB):
                        # Hash of the child, name of the child, parent:child relationship w/ parent ID in 1st position
        newAelements = [ ( md5(newA.encode('utf-8')).hexdigest(), newA, dataNewA[newA] ) for newA in dataNewA ]
        newBelements = [ ( md5(newB.encode('utf-8')).hexdigest(), newB, dataNewB[newB] ) for newB in dataNewB ]

        for hA, idA, dA in newAelements:
            for hB, idB, dB in newBelements:
                if hA < hB:                  
                    dX = dA
                    dY = dB
                else:            
                    dX = dB
                    dY = dA
                HoParameterSetObj = self.adjTree.getOrSet(idA, idB, HoParameterSet())
                #print(len(HoParameterSetObj))
                HoParameterSetObj.add(dX, dY)
                #print("#",len(HoParameterSetObj))

    def templateZipPair(self):
        templateColl = co.mdTree(append=False)
        for x,y,e in self.iterVisible():
            templates = e.templates
            #print(templates)
            for t1, t2 in zip(*templates):
                #print(t1, t2)
                templateColl.getOrSet(t1, t2, True)
        return templateColl

    def _weightProjector(self): # Classify element on which to project according to their number of 
                                # expected interactions in query network space
        #for domainElem in  self.hData:
        pass

    def T(self):
        pass

class HoParameterSet:
    def __init__(self):
        self.lowQueryParam = []
        self.highQueryParam = []
        self.visible = True
    def __repr__(self):
        return str({ "lowQueryParam" : [ p for p in self.lowQueryParam  if p.valid ],
                    "highQueryParam" : [ p for p in self.highQueryParam if p.valid ]
                })
    @property
    def depth(self):
        return len(self)
   
    @property
    def isEmpty(self):
        return len(self) == 0

    @property
    def templates(self):
        return ([ p.template for p in self.lowQueryParam  if p.valid ],
                [ p.template for p in self.highQueryParam  if p.valid ] 
                )
            
    def __len__(self):
        sz = 0
        for p in self.lowQueryParam:
            if p.valid:
                sz += 1
        return sz
   
    def add(self, x, y):
        self.lowQueryParam.append(  HoParameter(x) )
        self.highQueryParam.append( HoParameter(y) )
    
    def trim(self, simPct = 0.0, idPct = 0.0, cvPct = 0.0):
        self.visible = True
        for loHparam, hiHparam in self:
            #print (loHparam.simPct, loHparam.idPct, loHparam.cvPct)
            #print (hiHparam.simPct, hiHparam.idPct, hiHparam.cvPct)
            
            loHparam.valid = loHparam.simPct >= simPct and loHparam.idPct >= idPct and loHparam.cvPct >= cvPct
            hiHparam.valid = hiHparam.simPct >= simPct and hiHparam.idPct >= idPct and hiHparam.cvPct >= cvPct
            if not loHparam.valid or not hiHparam.valid:
                loHparam.valid = False
                hiHparam.valid = False
                #print("##OUT!!")

    def __iter__(self):
        for x,y in zip(self.lowQueryParam, self.highQueryParam):
            yield x,y
# Define properties for easy edge trimming
'''
()     [ 
      PSQ uniprot ID  
      R6_sequence Length,
      R6_hsp start position,
      R6_hsp stop position,
      PSQ_sequence Length,
      PSQ_hsp start position,
      PSQ_hsp stop position,
      HSP positive match number,
      HSP identical match number,
      HSP eValue
      ]
    '''
class HoParameter:
    def __init__(self, hVector):
        self.data = hVector
        self.valid = True
    
    def __repr__(self):
        return str(self.data)

    def __len__(self): #R6 aligned segment
        return int(self.data[3]) - int(self.data[2]) + 1
    
    @property
    def template(self):
        return self.data[0]
    @property
    def simPct(self):
        return  100.0 * float (self.data[7]) / float(len(self)) 
    @property
    def idPct(self):
        return 100.0 * float (self.data[8]) / float(len(self))
    @property
    def cvPct(self):
        return 100.0 * len(self) / int(self.data[1])
        
    @property
    def eValue(self):
        return self.data[8]
    




'''
A utility class to facilitate the visualisation of graph statistic in tabular format under jupyter
'''

class NodeView(object):
    def __init__(self, G):
        self.G = G
    def _repr_html_(self):
        htmlString = '<table><thead><th>Node</th><th>Degree</th></thead><tbody>'
        d = nx.degree(self.G)
        for n in sorted( [ { 'nodeRef' : k[0], 'degree' : k[1] } for k in d ], key=lambda x: x['degree'], reverse=True):
            color = 'brown' if self.G.node[n['nodeRef']]['group'] == 1 else 'black'

            htmlString += '<tr><td><a href="http://www.uniprot.org/uniprot/' + str(n['nodeRef']) + '" target="blank" style="color:' + color + '">' + str(n['nodeRef']) + '</td>'
            htmlString += '<td>' + str(n['degree']) + '</td></tr>'

        htmlString += '</tbody></table>'

        return htmlString