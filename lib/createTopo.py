import collections
import operator
import json
import re
#
# filePath
#
class Topology(object):

    def _sumPairs(self, dic):
        i = 0
        for k in dic:
            i += len(dic[k])
        return i

    @property
    def rawDicLen(self):
        return self._sumPairs(self.rawDic)

    @property
    def newDicLen(self):
        return self._sumPairs(self.newDic)

    def __init__(self):

        self.newDic = {}
        self.rawDic = {}
        self.sortedByOccurences = []
    allIdleft = set()
    orderedId = []
    keys = {}

    def merge_two_dicts(self, x, y):
        z = x.copy()
        z.update(y)
        return z

    def serialize(self, dictTopo, path):

        jsonStruct = {}
        for template, templateInInteraction in dictTopo.items():
            jsonStruct[template] = templateInInteraction

        json.dump(jsonStruct, file(path, 'w'))


    def deSerialize(self, path):

        with open (path, 'r') as file:
            data = json.load(file)

            # Renitialize the dictionnary if it's not already emtpty
            if not len(self.newDic) == 0:
                self.newDic.clear()

            for template, templateInInteraction in data.items():
                self.newDic[template] = templateInInteraction
        return self.newDic

    #Extract information from the database (Intact here)
    #Creat two dictionnary, one representing the full interaction database (dico)
    #The other one represent a reduced version of "dico", filtered with only interesting IDs
    #(ones which bring back QueryIDs with their BLAST)
    # First we keep pairs of uniprot identifier, discarding all other type of interactor pair
    def parseIntactMitab(self, filePath):
        dico = {}
        col1 = []
        col2 = []
        for line in open(filePath):
            sLine = line.split('\t')
            idOne = sLine[0]
            idTwo = sLine[1]
            idOne = idOne.split(':')
            idTwo = idTwo.split(':')
            if idOne[0] == "uniprotkb" and idTwo[0] == "uniprotkb":
                col1.append(idOne[1])
                col2.append(idTwo[1])

        CountIdOne = collections.Counter(col1)
        CountIdTwo = collections.Counter(col2)
        mergeIds = self.merge_two_dicts(CountIdOne, CountIdTwo)

        sortedByOccurencesDic = sorted(mergeIds.items(), key=operator.itemgetter(1), reverse=True)
        self.sortedByOccurences = [x[0] for x in sortedByOccurencesDic]

        colOne = list(col1)
        colTwo = list(col2)

        while len(colOne) > 0:
            for ids in self.sortedByOccurences:

                #print ids
                intercatWith = []
                toRemove = []
                i=0

                while i < len(colOne):

                    if colOne[i] == ids:
                        intercatWith.append(colTwo[i])
                        toRemove.append(i)
                    elif colTwo[i] == ids:
                        intercatWith.append(colOne[i])
                        toRemove.append(i)

                    i+=1

                    if intercatWith:
                        dico.update({ids : intercatWith})


                for toDel in sorted(toRemove, reverse=True):
                    del colOne[toDel]
                    del colTwo[toDel]

        self.rawDic = dico

    def filterWith(self, filterPath) :

        if(self.rawDicLen == 0):
            raise ValueError('You must parse a mitab file prior to filtering entries')

        allIdInR6 = []
        t = re.compile("^([\w]+)")
        for line in open(filterPath):
            if ":" not in line and line != "\n":
                m=t.match(line)
                if m:
                    allIdInR6.append(m.groups()[0])

        print("Number of intact uniprotID w/ R6 homologs " + str(len(allIdInR6)) )

        oldDic = self.rawDic.copy()
        tryId = ""
        for tryId in self.sortedByOccurences:
            if tryId not in allIdInR6:
                oldDic.pop(tryId, None)
                oldDic = {k: [e for e in v if e != tryId] for k, v in oldDic.iteritems()}

        self.newDic = dict((k, v) for k, v in oldDic.iteritems() if v)
        #return self.newDic