
class Divisome(object):
	def __init__(self):
		self.onlyId = []

	def getDivisomeID(self, path):
		with open (path, 'r') as file_div:
		    for i in xrange(1):
		        file_div.next()
		    for line in file_div:
		        sLine = line.split("\t")
		        self.onlyId.append(sLine[0])
		        
		return self.onlyId