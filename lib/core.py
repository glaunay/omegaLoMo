

import xml.etree.ElementTree as ET
import sys
import json
import numpy as np
from os import listdir
from os.path import isfile, join
import os
'''
.....
This implements the common element of the pipeline
....
'''



'''
 This describes the relationshionship between a query/template
'''
class Node(object):
    def __init__(self):
        self.template

    def __hash__(self):
        return hash(self.query)

    def __eq__(self, other):
        return self.query == other.query


class homologPair(object):
    def __init__(self, idTemplate, idQuery, param):
        self.template = idTemplate # A uniprot identifier or a unirpot object
        self.query = idQuery    # A uniprot identifier or a unirpot object
        self.param = param  ## store homology relationship quality descriptor %sim/%cov...

    def __str__(self):
        return self.query

    def serialize(self):
        pass

    def deSerialize(self):
        # read in a Json as string represrntaion of the object and populate its attributes
        pass

'''
This stores the colloction of homologs to a single Template element
'''
class hOmegaVector(object):
    def __init__(self, idTemplate= None, data = None):
        self.idTemplate = idTemplate
        self.data = data ## a list of homolPair objects

    def __hash__(self):
        return hash(self.idTemplate)

    def __str__(self):
        return self.idTemplate

    def _xmlRead(self, xmlFile, idQueryList):

        allQIdList = {}
        idList = []
        with open(idQueryList, 'r') as f:
            for line in f:
                idList.append(line.split("|")[1])
        fileName = xmlFile


        ## Extract from xml tree, the subtrees containing Hit_accession node text value present in idList
        ## Get last iteration
        tree=ET.parse(fileName)
        root = tree.getroot()
        parent_map = {c:p for p in root.iter() for c in p}

        blast_all_iter_node = root.find('./BlastOutput_iterations')
        self.idTemplate = root.find('./BlastOutput_query-def').text.split('|')[1]
        lastIter_subnode=root.findall('./BlastOutput_iterations/Iteration/Iteration_iter-num')[-1]
        lastIter_number=lastIter_subnode.text
        for iter_node in root.findall('./BlastOutput_iterations/Iteration'):
            if(iter_node.find('Iteration_iter-num').text != lastIter_number):
                #print "removing iteration " + str(iter_node)
                #print "from " + str(blast_all_iter_node)
                blast_all_iter_node.remove(iter_node)
                
        lastIter = parent_map[lastIter_subnode]
        Iteration_hits_node = lastIter.find("./Iteration_hits")
        for hit in Iteration_hits_node.findall("./Hit"):
            id = hit.find("Hit_accession")
            if id.text not in idList:
                #print "removing " + str(id) + ' from ' + str(Iteration_hits_node)
                Iteration_hits_node.remove(hit)
            else :
                allQIdList[id.text] = []
                coverList = []     
                ## Get hit score information and display stdout
                # loop over Hsp get Hsp_hit-from, Hsp_hit-to
                # permier arrive premier dedans
                # condition pour rentrer, etre non-chevauchant avec ceux deja presents.
                allowed = True
             

                #Iteration_hits_node.findall("./Hit/Hit_hsps/Hsp"):
                for hsp in hit.findall("./Hit_hsps/Hsp"):
                    Hfrom = hsp[6].text  
                    Hto = hsp[7].text
                    seq = ""

                    for seq in coverList:
                        #print seq, Hfrom, Hto
                  
                        if (int(seq[0]) <= int(Hfrom) <= int(seq[1])) or (int(seq[0]) <= int(Hto) <= int(seq[1])):
                            allowed = False
                        if (int(Hfrom) <= int(seq[0]) and  int(Hto) >= int(seq[1])):
                            allowed = False

                    if allowed:
                        allQIdList[id.text].append(([Hfrom, Hto, hit.find("Hit_hsps/Hsp/Hsp_positive").text, 
                                      hit.find("Hit_hsps/Hsp/Hsp_identity").text, 
                                      hit.find("Hit_hsps/Hsp/Hsp_evalue").text]))
                        
                        coverList.append((Hfrom, Hto))
        
        # Si on trouve des ID Query homologue a notre Template
        if len(allQIdList) > 0:
            for idQuery in allQIdList:
                if not self.data:
                    self.data = []
                self.data.append(homologPair(self.idTemplate, idQuery, allQIdList[idQuery]))
            return self.data
        else:
            return None 

    @property
    def empty(self): #True/False
        pass

    @property
    def id(self):
        return self.data[0].template

    def serialize(self):
        return { 'template ID' : self.idTemplate,
                 'Pairs Data' : [ d.__dict__ for d in self.data ]
                }

    def deSerialize(self):
        pass
'''
This is the collection of hOmegaVector
This has the dimension of the number of template with query hits
Guillaume would like us to implement an access by a query key !!


                


'''

#self.hOmegaData.add(xmlFile = blastPath, idQueryList=self.data['idQueryList'])


class HomegaSet(object):
    def __init__(self, **kwargs):
        self.data = [] # this sis the collection of hOmegaVector
        self.dict = {} # store array addres of item to be accessed through template id as a key

# Constructor mandatoray arguments
        if 'path' not in kwargs and 'bean' not in kwargs:
            raise ValueError ("Provide a blast results folders tree or a serialized omegaSet")
        if 'path' not in kwargs and 'queryIdList' not in kwargs:
             raise ValueError ("If you provide a blast results folders tree, you must supply a query list")
        else :
            self.idQueryList = kwargs["queryIdList"]

        if 'path' in kwargs:
            for xmlFileName in self._findXml(kwargs['path']):
                self.add(xmlFile=xmlFileName)
        
        elif 'bean' in kwargs:
            self.deSerialize(kwargs['bean'])
    
    def _findXml(self, path):
        if os.path.isdir(path):
            for rep in listdir(path):
                if not rep.startswith('.'):
                    blastPath = path + rep + "/logs/blast.out"
                    yield blastPath
        else:
            yield path

    def add(self, **kwargs):
        # item is a list or not , if not put it in a list
        # append item to
        if 'xmlFile' in kwargs:
            hOmegaVectorObj = hOmegaVector()
        if not hOmegaVectorObj._xmlRead( kwargs['xmlFile'], self.idQueryList):
            return

        if hOmegaVectorObj.idTemplate in self.dict:
            raise ValueError( str(hOmegaVectorObj.idTemplate) + " is already part of the set")
# Store current vetor in set
        self.data.append(hOmegaVectorObj)
        self.dict[hOmegaVectorObj.idTemplate] = hOmegaVectorObj
       
    def __getitem__(self, value):
        if value not in  self.dict:
            return None

        return self.dict[value]

    def OLD__getitem__(self, tup): ## return hOmegaVector object (or a list of ?? if query is the accessor)
        # x can be a sString or an index
        # A VOIR AVEC GUILLAUME (accesseur sous forme de matrix ou juste un Id de Template qui retourn les Querys)
        x, y = tup
        dataX = [hpair for hpair in self.dict[x].data]
        dataY = [hpair for hpair in self.dict[y].data]
        if len(dataX) > 0 :
            returnString = 'Template X : '+ dataX[0].template + ' <--> Queries X : ' + ' '.join([hpair.query for hpair in dataX])
        else :
            returnString += 'Empty Data for Template X'
        if len(dataY) > 0 :
            returnString += '\nTemplate Y : '+dataY[0].template + ' <--> Queries Y : ' + ' '.join([hpair.query for hpair in dataY])
        else :
            returnString += '\nEmpty Data for Template Y'
        return  returnString

    def __str__(self):
        return json.dumps({ "vectors" : [ v.serialize() for v in self.data ], 
                 "queryID" : self.idQueryList })

    def serialize(self, path):
        json.dump({ "vectors" : [ v.serialize() for v in self.data ], 
                 "queryID" : self.idQueryList }, file(path, 'w'))

    def deSerialize(self, path):
        with open (path, 'r') as file:
            data = json.load(file)
            fetch_data = [ hOmegaVector(v['template ID'], [homologPair( hp['template'], hp['query'], hp['param'])  for hp in v['Pairs Data']]) for v in data['vectors'] ]
            for vector in fetch_data:
                if vector.idTemplate in self.dict:
                    print 'Existe deja'
                else:
                    self.data.append(vector)
                    self.dict[vector.idTemplate] = vector

class fullMatrix(object):
    def __init__(self, data):
        ## self.data = { fullData: {}, idQueryList : [], blastFiles: [] }
        self.data = data
        self.hOmegaData = None

    def createSet(self):
        pass

    def reduce(self, HomegaSet):
        pass # return omegaMatrix

    def serialize(self):
        pass

    def deSerialize(self):
        pass

## implementaion would be a list of triplet of a pair of omegaVectors supplemned w/
# a free object containing structure or experimental information about the association between the two templates
##
class OmegaMatrix(object):
    def __init__(self, **kwargs):
        
        if not 'topo' in kwargs or not 'omegaSet' in kwargs:
            raise ValueError ('Provide a topographic file and a omegaSet')
        
        if 'topo' in kwargs and 'omegaSet' in kwargs:
            self.topo = kwargs['topo']
            self.omegaSet = kwargs['omegaSet']
            self.reduceTopo = {}
            self.dict = {}
            self.test = {}

    ## EXPLICATION
    def reduceAndVectorInject(self):
        
        if isinstance(self.topo, dict):
            for key, value in self.topo.iteritems():
                if key in self.omegaSet.dict:
                    self.reduceTopo[self.omegaSet.dict[key]] = []
                    self.test[key] = []

                    self.dict[key] = self.omegaSet.dict[key]
                else:
                    # Si la cle dans topo n'existe pas en tant que vecteur
                    continue

                for v in value:
                    if v in self.omegaSet.dict:
                        self.reduceTopo[self.omegaSet.dict[key]].append(self.omegaSet.dict[v])
                        self.test[key].append(v)

                        self.dict[v] = self.omegaSet.dict[v]
                    else:
                        continue

            #print self.dict
    
    def _multiMatrix(self, omegaVectorID_A, omegaVectorID_B):
        
        omegaVectorA = self.dict[omegaVectorID_A].data
        omegaVectorB = self.dict[omegaVectorID_B].data
        homologPairListA = [query for query in omegaVectorA]
        homologPairListB = [query for query in omegaVectorB]

        matrix = np.zeros((len(homologPairListA), len(homologPairListB)), dtype=object)

        for (x,y), value in np.ndenumerate(matrix):
            matrix[x, y] = [homologPairListA[x], homologPairListB[y]]

        #print str(homologPairListA) +'\n' +str(homologPairListB)+ '\n'+str(matrix)
        return matrix
    
    def project(self):
        pass

    def serialize(self):
        pass
    def deSerialize(self): ## Loic already has a serialized omega matrix as dict..
        pass

'''
    Each element of the queryMatrix
    (homologPair_I, homologPair_J, relationShipObject*)  # *:psciquic information relative to template pair OR Interface energy of original and threaded

'''
class queryMatrix(object):
    def __init__(self):
        pass

    def serialize(self):
        pass
    def deSerialize(self): ## Loic already has a serialized omega matrix as dict..
        pass

if __name__ == '__main__':
    data = {'repFile' : sys.argv[1],
       'idQueryList' : sys.argv[2]}

    omegaSet = HomegaSet(path=data['repFile'], queryIdList=data['idQueryList'])
    omegaSet.serialize(sys.argv[3])
