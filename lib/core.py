

import xml.etree.ElementTree as ET
import sys
'''
.....
This implements the common element of the pipeline
....
'''



'''
 This describes the relationshionship between a query/template
'''
class homologPair(object):
    def __init__(self, idTemplate, idQuery, param):
        self.template = idTemplate # A uniprot identifier or a unirpot object
        self.query = idQuery    # A uniprot identifier or a unirpot object
        self.param = param  ## store homology relationship quality descriptor %sim/%cov...

    def serialize(self):
        pass

    def deSerialize(self):
        # read in a Json as string represrntaion of the object and populate its attributes
        pass

'''
This stores the colloction of homologs to a single Template element
'''
class hOmegaVector(object):
    def __init__(self, idTemplate, blastFile, idQueryList):
        self.blastFile = blastFile
        self.idTemplate = idTemplate
        self.idQueryList = idQueryList
        self.data = [] ## a list of homolPair objects
        self.initialize()

    def initialize(self):
            #### Pouvoir recuperer Id de la Query
        resultParse = self._xmlRead(self.blastFile, self.idQueryList)
        for idQuery in resultParse.keys():
            self.data.append(homologPair(self.idTemplate, idQuery, resultParse[idQuery]))

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
        return allQIdList

    @property
    def empty(self): #True/False
        pass

    @property
    def id(self):
        return self.data[0].template

    def serialize(self):
        pass

    def deSerialize(self):
        pass
'''
This is the collection of hOmegaVector
This has the dimension of the number of template with query hits
Guillaume would like us to implement an access by a query key !!
'''

class hOmegaSet(object):
    def __init__(self, fullData):
        self.fulldata = fullData
        self.data = [] # this sis the collection of hOmegaVector
        self.dict = {} # store array addres of item to be accessed through template id as a key
        self.initialize()

    def initialize(self):
        # Pour chaque id template --> Vecteur --> Dict (Pour add seulement ?)
        for idQ in self.fulldata['blastFilesList'].keys():
            try:
                self.data.append(hOmegaVector(idQ , self.fulldata['blastFilesList'][idQ], self.fulldata['idQueryList'] ))
                self.dict[self.fulldata['fullData'].index(idQ)] = [idQ, self.data[-1]]
            except ValueError:
                print 'The vector cannot be created'
        print self.dict
    def add(self, **kwargs):
        # item is a list or not , if not put it in a list
        # append item to
        try:
            if not idTemplate in self.dict.keys():
                if 'xmlFile' in kwargs:
                    hOmegaVectorObj = hOmegaVector(xmlFile=kwargs['xmlFile'], idQueryList=self.idQueryList)
                    if not hOmegaVectorObj.empty:
                        self.data.append(hOmegaVectorObj)
                        # hOmegaVectorObj.id referencer dans le dict
        except ValueError:
            print 'We already know the queries id for the given template'

    def __getitem__(self, tup): ## return hOmegaVector object (or a list of ?? if query is the accessor)
        # x can be a sString or an index
        # A VOIR AVEC GUILLAUME (accesseur sous forme de matrix ou juste un Id de Template qui retourn les Querys)
        x, y = tup
        dataX = [hpair for hpair in self.dict[x][1].data]
        dataY = [hpair for hpair in self.dict[y][1].data]
        returnString = 'Template X : '+ dataX[0].template + ' <--> Queries X : ' + ' '.join([hpair.query for hpair in dataX])
        returnString += '\nTemplate Y : '+dataY[0].template + ' <--> Queries Y : ' + ' '.join([hpair.query for hpair in dataY])
        return  returnString

    def serialize(self):
        pass

    def deSerialize(self):
        pass

class fullMatrix(object):
    def __init__(self, data):
        ## self.data = { fullData: {}, idQueryList : [], blastFiles: [] }
        self.data = data
        self.hOmegaData = None

    def createSet(self):
        self.hOmegaData = hOmegaSet(self.data)
        return

    def reduce(self, hOmegaSet):
        pass # return omegaMatrix

    def serialize(self):
        pass

    def deSerialize(self):
        pass

## implementaion would be a list of triplet of a pair of omegaVectors supplemned w/
# a free object containing structure or experimental information about the association between the two templates
##
class omegaMatrix(object):
    def __init__(self):
        pass
    #self.data=[ (Hvec, Hvec, relationShipObject), ... () ]
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
