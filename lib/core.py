import xml.etree.ElementTree as ET
import sys
import json
import numpy as np
from os import listdir
from os.path import isfile, join
import os
import fnmatch
import gzip
'''
.....
This implements the common element of the pipeline
....
'''



'''
 This describes the relationshionship between a query/template
'''

''' This class represent the interaction between a template and a query and can be used by networkX as a node'''
class Node(object):
    def __init__(self, homologPair):
        self.template = homologPair.template
        self.query = homologPair.query
        self.param = homologPair.param

    def __hash__(self):
        return hash(self.query)

    def __str__(self):
        return self.query

    def __eq__(self, other):
        return self.query == other.query

    def __gt__(self, other):
        if cmp(self.query, other.query) > 0 :
            return True
        return False

    def __lt__(self, other):
        if cmp(self.query, other.query) < 0:
            return True
        return False
    def __repr__(self):
        return self.query


''' This class store a node object which define the homology between a template and a query'''
class homologPair(object):
    def __init__(self, idTemplate, idQuery, param):
        self.template = idTemplate # A uniprot identifier or a unirpot object
        self.query = idQuery    # A uniprot identifier or a unirpot object
        self.param = param  ## store homology relationship quality descriptor %sim/%cov...

    # The string representation of homologPair object return th query Id
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

    def __len__(self):
        return len(self.data)

    def __hash__(self):
        return hash(self.idTemplate)

    def __str__(self):
        return self.idTemplate

    # This method take as input a XML Blast Output (xmlFile) that can filter using a query list (idQueryList)
    def _xmlRead(self, xmlFile, idQueryList):

        getXML = []
        currentXML = None

        f_reader = open
        if xmlFile.endswith("gz"):
            f_reader = gzip.open


        for item in self.xmlSplitter(f_reader(xmlFile)):
            doc = ET.fromstring(item)
            length = doc.find('./BlastOutput_query-len').text
            getXML.append((length, doc))
        currentXML = [xml for xml in getXML if xml[0] == max(length[0] for length in getXML if length != 0)]

        allQIdList = {}
        idList = []

        with open(idQueryList, 'r') as f:
            for line in f:
                idList.append(line.split("|")[1])
        fileName = xmlFile

        ## Extract from xml tree, the subtrees containing Hit_accession node text value present in idList
        ## Get last iteration
        #tree=ET.parse(fileName)
        #tree.getroot()
        root = currentXML[0][1]
        parent_map = {c:p for p in root.iter() for c in p}

        blast_all_iter_node = root.find('./BlastOutput_iterations')
        self.idTemplate = root.find('./BlastOutput_query-def').text.split('|')[1]
        lenQuery = root.find('./BlastOutput_query-len').text
        lastIter_subnode=root.findall('./BlastOutput_iterations/Iteration/Iteration_iter-num')[-1]
        lastIter_number=lastIter_subnode.text
        for iter_node in root.findall('./BlastOutput_iterations/Iteration'):
            if(iter_node.find('Iteration_iter-num').text != lastIter_number):
                #print "removing iteration " + str(iter_node)
                #print "from " + str(blast_all_iter_node)
                blast_all_iter_node.remove(iter_node)

        lastIter = parent_map[lastIter_subnode]
        Iteration_hits_node = lastIter.find("./Iteration_hits")
        if Iteration_hits_node is not None:
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
                            allQIdList[id.text].append(([Hfrom,
                                          Hto,
                                          hit.find("Hit_hsps/Hsp/Hsp_positive").text,
                                          hit.find("Hit_hsps/Hsp/Hsp_identity").text,
                                          hit.find("Hit_hsps/Hsp/Hsp_evalue").text,
                                          lenQuery,
                                          hit.find("Hit_hsps/Hsp/Hsp_query-from").text,
                                          hit.find("Hit_hsps/Hsp/Hsp_query-to").text,
                                          hit.find("Hit_len").text]))

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

    def printHomologs(self):
        for homolPairObj in self.data:
            #print homolPairObj.param
            print homolPairObj.query + ":" + '|'.join(','.join(s) for s in homolPairObj.param) + '\n'

    def xmlSplitter(self, data, separator=lambda x: x.startswith('<?xml')):
        buff = []
        for line in data:
            if separator(line):
                if buff:
                    yield ''.join(buff)
                    buff[:] = []
            buff.append(line)
        yield ''.join(buff)


    @property
    def id(self):
        return self.data[0].template

    # Generate a dictionnary format which sums up all homologPair object of the current vector (JSon format)
    def serialize(self):
        return { 'template ID' : self.idTemplate,
                 'Pairs Data' : [ d.__dict__ for d in self.data ]
                }

    def deSerialize(self):
        pass


'''
This is the collection of hOmegaVector
This has the dimension of the number of template with query hits
'''
class HomegaSet(object):
    def __init__(self, **kwargs):
        self.data = [] # this sis the collection of hOmegaVector
        self.dict = {} # store array addres of item to be accessed through template id as a key
        self.idQueryList = None

# Constructor mandatoray arguments
        if 'queryIdList' in kwargs:
            self.idQueryList = kwargs["queryIdList"]

        if 'blastXmlPath' in kwargs:
            print("path is here " + kwargs['blastXmlPath'])
            for xmlFileName in self._findXml(kwargs['blastXmlPath']):
                print("reading from " + xmlFileName)
                self.add(xmlFile=xmlFileName)

        if 'jsonSerialPath' in kwargs:
            print("path is here " + kwargs['jsonSerialPath'])
           # jsonList = [ fJson for fJson in self._findJson(kwargs['jsonSerialPath']) ]

            self.deSerialize( jsonIndividualFileList = self._findJson(kwargs['jsonSerialPath']) )
                #self.add(jsonFile=jsonFileName)


        elif 'bean' in kwargs:
            self.deSerialize(beanPath=kwargs['bean'])

    # Check if the path is a repository of JSon file or single file
    # Browse recursively the path if is a repository and catch all JSon files

    # Check if the path is a repository of XML file or single file
    # Browse recursively the path if is a repository and catch all XML files
    def _findXml(self, path):
        if os.path.isdir(path):
            for root, subdirs, filenames in os.walk(path):
                for filename in fnmatch.filter(filenames, '*.out*'):
                    yield os.path.join(root, filename)
        else:
            if path.endswith('.out') or path.endswith('.out.gz'):
                yield path
    def _findJson(self, path):
        if os.path.isdir(path):
            for root, subdirs, filenames in os.walk(path):
                for filename in fnmatch.filter(filenames, '*.json*'):
                    yield os.path.join(root, filename)
        else:
            if path.endswith('.json') or path.endswith('.json.gz'):
                yield path

    # This method can take as arguments an XML Blast Output
    # Check if the hOmegaVector generated by the previous file already exist in the omegaSet
    # If not exist It will be append in th current omegaSet
    def add(self, **kwargs):
        hOmegaVectorObj = hOmegaVector()

        # item is a list or not , if not put it in a list
        # append item to
        if 'xmlFile' in kwargs:
            if not hOmegaVectorObj._xmlRead( kwargs['xmlFile'], self.idQueryList):
                return
        if 'jsonFile' in kwargs:
            if not hOmegaVectorObj._jsonRead( kwargs['jsonFile'], self.idQueryList):
                return


        if hOmegaVectorObj.idTemplate in self.dict:
            raise ValueError( str(hOmegaVectorObj.idTemplate) + " is already part of the set")

        # Store current vector in set
        self.data.append(hOmegaVectorObj)
        self.dict[hOmegaVectorObj.idTemplate] = hOmegaVectorObj

    def __getitem__(self, value):
        if value not in  self.dict:
            return None

        return self.dict[value]

    def __str__(self):
        jsonStruct = { "vectors" : [ v.serialize() for v in self.data ], "queryID" : None}
        if 'self.idQueryList' in globals():
            jsonStruct[queryID] = self.idQueryList

        return json.dumps(jsonStruct)

    def __iadd__(self, other):
        for vector in other.data:
            if not vector.idTemplate in self.dict:
                self.data.append(vector)
                self.dict[vector.idTemplate] = hOmegaVectorObj
                return self
            else:
                raise ValueError ("The vector is "+vector.idTemplate+" already in the omegaSet")

    def __add__(self, other):
        for vector in other.data:
            if not vector.idTemplate in self.dict:
                self.data.append(vector)
                self.dict[vector.idTemplate] = hOmegaVectorObj
                return self
            else:
                raise ValueError ("The vector is "+vector.idTemplate+" already in the omegaSet")

    def __len__(self):
        return len(self.data)

    # Generated a JSon file
    def serialize(self, path):
        jsonStruct = { "vectors" : [ v.serialize() for v in self.data ], "queryID" : None}
        if self.idQueryList:
            jsonStruct['queryID'] = self.idQueryList

        json.dump(jsonStruct, file(path, 'w'))


    def deSerialize(self, beanPath=None, jsonIndividualFileList = None):
        # We go through all vectors element w/ a single dictionary/serialized dataset and yield each one
        def _iterAll (data):
            for v in data['vectors']:
                yield v
        # We go through all specified json files lookking for vector element and yiel each one vectors element w/ a single dictionary/serialized dataset

        def _iterOnJsonFile(fileListIter):
            for fName in fileListIter:
        #        print("Opening " + fName)
                with open (fName, 'r') as file:
                    data = json.load(file)
                    for v in data['vectors']:
                        yield v

        inputStuff = None
        iter_anon = None
        if beanPath :
            with open (beanPath, 'r') as file:
                inputStuff = json.load(file)
            iter_anon = _iterAll

        elif jsonIndividualFileList:
            iter_anon = _iterOnJsonFile
            inputStuff = jsonIndividualFileList
        else:
            raise TypeError("Full Bean of list of json files required for deserialization")
            #for v in data['vectors']:
        for v in iter_anon(inputStuff):
                #print v
            if not v['template ID'] in self.dict:
               # print "Adding new template ID at " + v['template ID']
                vector = hOmegaVector(v['template ID'])
                    #print vector
                if vector.data is None:
                    vector.data = []
                for hp in v['Pairs Data']:
                    vector.data.append(homologPair( hp['template'], hp['query'], hp['param']))
                self.data.append(vector)
                self.dict[vector.idTemplate] = vector

            else:
                print (v['template ID'] + ' Already exist')
               # print ("-->" + str(self.dict[v['template ID']]))
               # return;


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
            self.templatePairs = []
            self.queryTopo = {} ### ATTRIBUT ?

    def reduceAndVectorInject(self):

        if isinstance(self.topo, dict):
            for key, value in self.topo.iteritems():
                if key in self.omegaSet.dict:
                    self.reduceTopo[self.omegaSet.dict[key]] = []
                    self.dict[key] = self.omegaSet.dict[key]
                else:
                    # Si la cle dans topo n'existe pas en tant que vecteur
                    continue

                for v in value:
                    if v in self.omegaSet.dict:
                        self.reduceTopo[self.omegaSet.dict[key]].append(self.omegaSet.dict[v])
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
            matrix[x, y] = (Node(homologPairListA[x]), Node(homologPairListB[y]))

        #print str(homologPairListA) +'\n' +str(homologPairListB)+ '\n'+str(matrix)
        return matrix

    def project(self): ## C'est project qui doit prenre un fichier de liste de pair ou _getTemplatePairs ?
        queryMatrixObj = QueryMatrix()
        for (omegaVector_A, omegaVector_B) in self._getTemplatePairs():
            queryTopology = self.mapMiniMatrix(omegaVector_A, omegaVector_B)
            queryMatrixObj.add(queryTopology)

        return queryMatrixObj

# Generate all pairs of template objects in interaction

    def _getTemplatePairs(self):
        for pKey in self.reduceTopo:
            for sKey in self.reduceTopo[pKey]:
                yield (pKey, sKey)

    def mapMiniMatrix(self, omegaVector_A, omegaVector_B):
        miniMatrix = self._multiMatrix(omegaVector_A.idTemplate, omegaVector_B.idTemplate)
        queryTopo = {}
        for x in range(miniMatrix.shape[0]):
            for y in range(miniMatrix.shape[1]):
                (NodeObjA, NodeObjB) = miniMatrix[x, y]

                # These 2 queries in interactions
                #queryID_A = pos[0]
                #queryID_B = pos[1]

                #storeArray = None

                lo_query = NodeObjA if NodeObjA < NodeObjB else NodeObjB # Opt-in 1ry key
                hi_query = NodeObjA if NodeObjA > NodeObjB else NodeObjB
                if NodeObjA  == NodeObjB:
                    lo_query = NodeObjA
                    hi_query = NodeObjB

                if lo_query in queryTopo:
                    #print queryTopo[lo_query]

                    if hi_query not in queryTopo[lo_query]:
                        queryTopo[lo_query][hi_query] = []
                else:
                    queryTopo[lo_query] = { hi_query : [] }


                queryTopo[lo_query][hi_query].append({'loQueryParam' : lo_query.param[0], 'hiQueryParam' : hi_query.param[0]})

                #DEPRECATED
                #storeArray = queryTopo[lo_query][hi_query]
                # FORMAT PARAMETERS ?
                #storeArray.append({'loQueryParam' : lo_query.param[0], 'hiQueryParam' : hi_query.param[0]})
            #print queryTopo
        return queryTopo

    def serialize(self):
        pass
    def deSerialize(self): ## Loic already has a serialized omega matrix as dict..
        pass


class QueryMatrix(object):
    def __init__(self):
        self.queryTopo = []
        self.ghost_edges = []
        self.dictQuery = {}

    # Refreshing global topology with a subtopoloy, aka miniMatrix obtained by cartesian product of 2 omegaVectors
    def add(self, queryTopo):
        self.queryTopo.append(queryTopo)

        # Add memories adresses for a given Node
        for nKey, nVal in queryTopo.iteritems():
           # print str(nKey) + " ___ " + str(nVal)
            if not nKey.query in self.dictQuery:
                self.dictQuery[nKey.query] = [nKey]
            else:
                self.dictQuery[nKey.query].append(nKey)

            for hQuery in nVal:
                if not hQuery.query in self.dictQuery:
                    self.dictQuery[hQuery.query] = [hQuery]
                else:
                    self.dictQuery[hQuery.query].append(hQuery)


    def getEdges(self, **kwargs):

        for interaction in self.queryTopo:
            for lowQuery in interaction:
                for highQuery in interaction[lowQuery]:
                    yield {'lowQuery' : lowQuery, 'highQuery' : highQuery,
                    'loQueryEval' : [param['loQueryParam'] for param in interaction[lowQuery][highQuery]],
                    'hiQueryEval' : [param['hiQueryParam'] for param in interaction[lowQuery][highQuery]]}

    def serialize(self):
        pass

    def deSerialize(self): ## Loic already has a serialized omega matrix as dict..
        pass

# For sbatch blast post-processing


# GL -- 12/20/2017
# Stand-alone call of the form
# python core.py --blast2log/--blast2json  BLAST_FILE ID_QUERY_LIST_FILE
if __name__ == '__main__':
    switch = sys.argv[1]
    data = {'repFile' : sys.argv[2],
       'queryIdList' : sys.argv[3]}

    if switch == '--blast2log':
        hOmegaVectorObj = hOmegaVector()
        hOmegaVectorObj._xmlRead(data['repFile'], data['queryIdList'])
        hOmegaVectorObj.printHomologs()

    elif switch == '--blast2json':
        omegaSet = HomegaSet(blastXmlPath=data['repFile'], queryIdList=data['queryIdList'])
        omegaSet.serialize(sys.argv[4])
    else :
        raise ValueError('unrecognized command switch' + switch)