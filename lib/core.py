'''
.....
This implements the common element of the pipeline
....
'''



'''
 This describes the relationshionship between a query/template
'''
class homologPair(object):
   def __init__(self):
        self.template = None # A uniprot identifier or a unirpot object
        self.query = None    # A uniprot identifier or a unirpot object
        self.param = {} ## store homology relationship quality descriptor %sim/%cov...
        pass

    def serialize(self):
        # return a Json as string represrntaion of the object
        pass

    def deSerialize(self):
        # read in a Json as string represrntaion of the object and populate its attributes
        pass

'''
This stores the colloction of homologs to a single Template element
'''
class hOmegaVector(object):
    def __init__(self):
        self.data = [] ## a list of homolPair objects
        pass

    def _xmlRead(xmlFile, idQueryList):
        pass ## xmlParser

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
    def __init__(self, idQueryList):
        self.idQueryList = idQueryList
        self.data = []# this sis the collection of hOmegaVector
        self.dict = {} # store array addres of item to be accessed through template id as a key
    def add(self, **kwargs):
        # item is a list or not , if not put it in a list
        # append item to
        if 'xmlFile' in kwargs:
            hOmegaVectorObj = hOmegaVector(xmlFile=kwargs['xmlFile'], idQueryList=self.idQueryList)
            if !hOmegaVectorObj.empty:
                self.data.append(hOmegaVectorObj)
                # hOmegaVectorObj.id referencer dans le dict
                # raise error if already present
    def __getitem__(x): ## return hOmegaVector object (or a list of ?? if query is the accessor)
        # x can be a sString or an index
        return self.data[x]

    def serialize(self):
        pass

    def deSerialize(self):
        pass

class fullMatrix(object):
    __init__(self):
        pass

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
    __init__(self):
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
