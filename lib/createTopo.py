import collections
import operator
import json


class Topology(object):
    def __init__(self, filePath, filterPath):
        self.filePath = filePath
        self.filterPath = filterPath
        self.newDic = {}
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
        for template, templateInInteraction in dictTopo.iteritems():
            jsonStruct[template] = templateInInteraction
            
        json.dump(jsonStruct, file(path, 'w'))


    def deSerialize(self, path):
    
        with open (path, 'r') as file:
            data = json.load(file)
            
            # Renitialize the dictionnary if it's not already emtpty
            if not len(self.newDic) == 0:
                self.newDic.clear()
            
            for template, templateInInteraction in data.iteritems():
                self.newDic[template] = templateInInteraction
        return self.newDic

    #Extract information from the database (Intact here)
    #Creat two dictionnary, one representing the full interaction database (dico)
    #The other one represent a reduced version of "dico", filtered with only interesting IDs 
    #(ones which bring back QueryIDs with their BLAST)
    def parse_Intact_mitab(self, filePath):
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
                    
        return dico

    def filter_With(self, filterPath) :
        
        dico = self.parse_Intact_mitab(self.filePath)
        
        allIdInR6 = []
        for line in open(filterPath):
            if ":" not in line and line != "\n":
                allIdInR6.append(line[:-1])

        oldDic = dico.copy()
        tryId = ""
        for tryId in self.sortedByOccurences:
            if tryId not in allIdInR6:
                oldDic.pop(tryId, None)
                oldDic = {k: [e for e in v if e != tryId] for k, v in oldDic.iteritems()}

        self.newDic = dict((k, v) for k, v in oldDic.iteritems() if v)
        return self.newDic