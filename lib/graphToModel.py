import urllib3
import os
import fnmatch
import re
import json

import pyproteins.homology.dummyThreaders as thr
import pyproteinsExt.structure.coordinates as PDB
import pyproteins.homology.query as qr
import pyproteins.services.utils

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

	'''
			CLASS : GRAPH TO MODEL

		INPUT:
			None

		ARGUMENTS:
			- graphObj : The graph object used to get all possible models
			- modelForEdge : Dictionnary where all models for a edge is refered
			- _modelGenerated : Refered the generate models for the edge object 
			- nodeList : All nodes used to generate models provide by the user
	'''

class GraphToModel(object):
	def __init__(self):
		self.graphObj = None
		self.modelForEdge = {}
		self._modelGenerated = {}
		self.nodeList = []

	'''
	
	HOMOLOGUE IN BLAST

	Input : 
	- edgeObj : The given edge object
	- Node : Node object
	- pathBlast : Path to folders containing blast files
	- pdb_files : Path to PDB files of monomeric chains of InterEvol DB
	
	Output : None

	Function : 

	edgeObj and Node can't be filled in the same time. If edgeObj is provided
	the two nodes will get seperated in a "right" vertice and "left" vertice. Each one
	corresponding to a blast output file. For a blast file, the function will fetch all
	the hits with :
	- the monomeric ID of the structural dimer
	- the monomeric chain of the strutural dimer
	- the e-value corresponding to the hit. 
	These hits will be stock in the edge Object in two "Hits" list corresponding
	to the "right" vertice and the "left" vertice.

	If Node provide instead of edgeObj, this will do the same job tha previously. It will
	fetch all hits with information describe above in the blast output corresponding to the node file and stock these in
	self.node argument of the current class.

	'''

	def homologuesInBlast(self, edgeObj, Node, pathBlast , pdb_files):

		if edgeObj != None and Node == None:

			verticeR = edgeObj.v1.query
			verticeL = edgeObj.v2.query

			verticePathBlast = pathBlast+self.graphObj.queryCenter+"/"

			verticeRBlast = [blast for blast in self._findBlast(verticeR, verticePathBlast)][0]
			verticeLBlast = [blast for blast in self._findBlast(verticeL, verticePathBlast)][0]
			
			filesVerticeR = verticeRBlast.split(".")[0]
			filesVerticeL = verticeLBlast.split(".")[0]

			for verticeRitem in self._blastSplitter(open(verticeRBlast)):
				#fromNode = blast.split("/")[5]

				verticeR_Evalue = re.search('Expect =[ .0-9a-zA-Z_-]+', verticeRitem)
				verticeR_idChain = re.search('^>[A-Za-z0-9_]+', verticeRitem)
				
				if verticeR_Evalue:
					verticeR_Evalue = verticeR_Evalue.group(0).replace(" ", "")
					
					if verticeR_Evalue.split("=")[1].startswith("e"):
						verticeR_Evalue = float("1"+verticeR_Evalue.split("=")[1])
					else:
						verticeR_Evalue = float(verticeR_Evalue.split("=")[1])
				else:
					pass
					#Maybe no hits in .BLAST output

				if verticeR_idChain:
					verticeR_Struct = verticeR_idChain.group(0).split(">")[1]
					verticeR_idChain = verticeR_Struct.split("_")[0]
					verticeR_Chain = verticeR_Struct.split("_")[1]
				else:
					pass
					#Maybe no hits in .BLAST output
				if verticeR_Evalue and verticeR_idChain:
					edgeObj.hitsV1.append((verticeR, verticeR_idChain, verticeR_Chain, verticeR_Evalue, filesVerticeR, pdb_files))

			for verticeLitem in self._blastSplitter(open(verticeLBlast)):

				verticeL_Evalue = re.search('Expect =[ .0-9a-zA-Z_-]+', verticeLitem)
				verticeL_idChain = re.search('^>[A-Za-z0-9_]+', verticeLitem)
				
				if verticeL_Evalue:
					verticeL_Evalue = verticeL_Evalue.group(0).replace(" ", "")

					if verticeL_Evalue.split("=")[1].startswith("e"):
						verticeL_Evalue = float("1"+verticeL_Evalue.split("=")[1])
					else:
						verticeL_Evalue = float(verticeL_Evalue.split("=")[1])
				else:
					pass
					#Maybe no hits in .BLAST output

				if verticeL_idChain:
					verticeL_Struct = verticeL_idChain.group(0).split(">")[1]
					verticeL_idChain = verticeL_Struct.split("_")[0]
					verticeL_Chain = verticeL_Struct.split("_")[1]
				else:
					pass
					#Maybe no hits in .BLAST output

				if verticeL_Evalue and verticeL_idChain:
					edgeObj.hitsV2.append((verticeL, verticeL_idChain, verticeL_Chain, verticeL_Evalue, filesVerticeL, pdb_files))
		
		############# NODE PART ##############
		else:
			node = Node.query
			nodePathBlast = pathBlast+node+"/"

			nodeBlast = [blast for blast in self._findBlast(node, nodePathBlast)][0]
			filesNode = nodeBlast.split(".")[0]

			for nodeItem in self._blastSplitter(open(nodeBlast)):
				node_Evalue = re.search('Expect =[ .0-9a-zA-Z_-]+', nodeItem)
				node_idChain = re.search('^>[A-Za-z0-9_]+', nodeItem)
				
				if node_Evalue:
					node_Evalue = node_Evalue.group(0).replace(" ", "")
					
					if node_Evalue.split("=")[1].startswith("e"):
						node_Evalue = float("1"+node_Evalue.split("=")[1])
					else:
						node_Evalue = float(node_Evalue.split("=")[1])
				else:
					pass
					#Maybe no hits in .BLAST output

				if node_idChain:
					node_Struct = node_idChain.group(0).split(">")[1]
					node_idChain = node_Struct.split("_")[0]
					node_Chain = node_Struct.split("_")[1]
				else:
					pass
					#Maybe no hits in .BLAST output
				if node_Evalue and node_idChain:
					Node.monomer.append((node, node_idChain, node_Chain, node_Evalue, filesNode, pdb_files))

			self.nodeList.append(Node)
				

	'''
	
	CREATE SET OF DIMERRIC

	Input : 
	- edgeObj : The Edge object containing two nodes
	- eValue : Float e-value treshold  
	
	Output :
	Raise a tuple of informations:
	- id_Dimeric : Id of the dimeric reference structure
	- chain_MonomerA , chain_MonomerB : The monomeric chain corresponding to the monomers invested in a demeric structure
	- pathFiles_DimericA, pathFiles_DimericB : Path to the monomeric PDB file invested in the same dimeric structure
	- pdb_files : Path to all PDB files of dimeric reference structures
	- evalue_MonomerA, evalue_MonomerB : The e-vale corresponding to the monomeric hits invested in the same dimeric reference structure

	Function:
	This method find all hits for vertices invested in an edge under the e-value.
	It generat a set of possible combination between the monomeric hits of the two vertices
	'''

	def _createSetOfDimeric(self, edgeObj, eValue):
		setPotentialDimeric = set()
			
		for structA in edgeObj.hitsV1:
			if float(structA[3]) <= float(eValue):
				id_MonomerA = structA[1]
				chain_MonomerA = structA[2]

			for structB in edgeObj.hitsV2:
				if float(structB[3]) <= float(eValue):
					id_MonomerB = structB[1]
					chain_MonomerB = structB[2]

					#struct[0] --> id Query 
					if structA[0] != structB[0]:
						if id_MonomerA == id_MonomerB:
							id_Dimeric = id_MonomerA
							evalue_MonomerA = structA[3]
							evalue_MonomerB = structB[3]
							pathFiles_DimericA = structA[4]
							pathFiles_DimericB = structB[4]
							pdb_files = structA[5]

							yield(id_Dimeric, chain_MonomerA, chain_MonomerB, pathFiles_DimericA, pathFiles_DimericB, pdb_files, evalue_MonomerA, evalue_MonomerB)

	'''
	
	GENERATE TEMPLATE

	Input :
	- fullTopo : A topography pf a full dataset like a proteom
	- fullGraph : The graph generated from the fullTopo
	- listCentralNodes : List of nodes used to generate a graph of the first neighbors 
	- dataBin : Dictionnary containing the necessary path to the binaries
	- eValue : Float e-value

	Output : None

	Function:
	--> Method call in first
	For each nodes in listCentralNodes It will create the first neighbors graph from the fullGraph.
	For each neighboring node the method will execute the self.homologuesInBlast() and self.findDimericTemplates
	'''

	def generateTemplates(self, fullTopo, fullGrah, listCentralNodes, dataBin, eValue):

		# Verification step on binarie files
		if not "hhbin" in dataBin:
			raise valueError("Provide a path to the HHSuite Binaries, e.g : ~/Desktop/pipeline/hhsuite/bin/ ")
		elif not "workDir" in dataBin:
			raise valueError("Provide a filepath which will contain the models, e.g : ~/Desktop/Modelisation/")
		elif not "pathBlast" in dataBin:
			raise valueError("Provide a filepath which will contain the blast files, e.g : ~/Desktop/fastaGraph/")
		elif not "pdb_files" in dataBin:
			raise valueError("Provide a filepath which will contain the pdb files, e.g : ~/Desktop/interevol/InterEvol/PDB/interEvol_structures")
		elif not "path_table_dimeric_struct" in dataBin:
			raise valueError("Provide a filepath which will contain the dimeric association chains, e.g : ~/Desktop/interevol/InterEvol/PDB/interEvol_structures")

		#Create a neighbor nodes graph for each protein in listCentralNodes
		for centralNode in listCentralNodes:

			premierVoisin = fullTopo.createNeiGraph(centralNode, fullGrah)
			generateNeighboor = premierVoisin.generateNeighboorGraph
			self.graphObj = premierVoisin

		for nodeR in premierVoisin.graphData.edge:
			if nodeR.query == centralNode:
				for nodeL in premierVoisin.graphData.edge[nodeR]:
					edge = Edge(nodeR, nodeL)
					self.homologuesInBlast(edge, None, dataBin["pathBlast"], dataBin["pdb_files"])
					self.findDimericTemplates(edge, 1.0, dataBin["hhbin"], dataBin["workDir"], dataBin["path_table_dimeric_struct"])

	'''
	
	FIND DIMERIC TEMPLATES

	Input : 
	- edgeObj : The edge object
	- eValue : float e-value
	- hhbin : The path to the HHSuite binaries
	- workDir : The path where the data will be generated
	- path_table_dimeric_struct : The path to the file used to check all monomeric associations

	Output : 
	- None

	Function : 
	Get all combinations generated bu the self._createSetOfDimeric method.
	For a combination of monomeric strcures It will check if this one is refered if the table (path_table_dimeric_struct) as a dimeric structure.
	If It's refered this lunch a threading operation.

	It calls also self._graphModelPerVertice method 
	'''


	def findDimericTemplates(self, edgeObj, eValue, hhbin, workDir, path_table_dimeric_struct):
		
		setPotentialDimericStruct = set([dimeric for dimeric in self._createSetOfDimeric(edgeObj, eValue)])
		workDir += self.graphObj.queryCenter+'/'+edgeObj.v1.query+'_'+edgeObj.v2.query+'/'
		model = 0
		evalue_Dimeric = {}

		with open(path_table_dimeric_struct, 'r') as table:
			for line in table:
				for dimericStructure in setPotentialDimericStruct:
					
					id_struct = dimericStructure[0]
					chaine_structA = dimericStructure[1]
					chaine_structB = dimericStructure[2]
					pathFiles_structA = dimericStructure[3]
					pathFiles_structB = dimericStructure[4]
					pdb_files = dimericStructure[5]
					evalue_MonomerA = dimericStructure[6]
					evalue_MonomerB = dimericStructure[7]

					if line.startswith(id_struct+"\t"+chaine_structA+"-"+chaine_structB):

						# If dimeric structure exist in PDB format
						if os.path.isfile(pdb_files+"/"+id_struct+"_"+chaine_structA+".pdb") and os.path.isfile(pdb_files+"/"+id_struct+"_"+chaine_structB+".pdb"):

							if not edgeObj.v1.query in self._modelGenerated:
								self._modelGenerated[edgeObj.v1.query] = {edgeObj.v2.query : { }}
							
							if not edgeObj.v2.query in self._modelGenerated[edgeObj.v1.query]:
								self._modelGenerated[edgeObj.v1.query][edgeObj.v2.query] = { }

							if not id_struct in self._modelGenerated[edgeObj.v1.query][edgeObj.v2.query]:
								self._modelGenerated[edgeObj.v1.query][edgeObj.v2.query][id_struct] = {
									"pir" : workDir+"/Default.pir",
									"pdb_monom_R" : pdb_files+"/"+id_struct+"_"+chaine_structA+".pdb",
									"pdb_monom_L" : pdb_files+"/"+id_struct+"_"+chaine_structB+".pdb",
									"query_monom_R" : pathFiles_structA+".fasta",
									"query_monom_L" : pathFiles_structB+".fasta",
									"evalue_monom_R" : str(evalue_MonomerA),
									"evalue_monom_L" : str(evalue_MonomerB)
									}

								# Add id Struct folder
								workDirFinal = workDir
								workDirFinal += id_struct

								self._launchModelisation(edgeObj, id_struct, [chaine_structA, chaine_structB], [evalue_MonomerA, evalue_MonomerB], [pathFiles_structA+".fasta", pathFiles_structB+".fasta"], [pdb_files+"/"+id_struct+"_"+chaine_structA+".pdb", pdb_files+"/"+id_struct+"_"+chaine_structB+".pdb"], workDirFinal, hhbin)

								evalue_Dimeric[id_struct] = [(id_struct+"_"+chaine_structA, evalue_MonomerA), 
								(id_struct+"_"+chaine_structB, evalue_MonomerB)]
								
								model += 1
			
			if edgeObj.v1.query not in self.modelForEdge:
				self.modelForEdge[edgeObj.v1.query] = {edgeObj.v2.query : [model, evalue_Dimeric]}
			else:
				self.modelForEdge[edgeObj.v1.query][edgeObj.v2.query] = [model, evalue_Dimeric]
			
			if len(self._modelGenerated) > 0:
				self._graphModelPerVertice(edgeObj)
				self._modelGenerated.clear()

	'''
	
		SERIALIZE

	INPUT : 
	- data : Data as dictionnary used to serialize
	- path : The path where the JSON file will be generate

	OUTPUT: None

	FUNCTION : 
	Create a JSON file by writing the data inside

	'''

	def _serialize(self, data, path):
		print "Serialiaze"

		json.dump(data, file(path, 'w'))

	'''
		LUNCH MODELISATION

	INPUT :
		- edgeobj : The edge Object
		- id_struct : The ID of the PDB dimeric structure  
		- chains : List of monomerics chains involved in the same dimeric structure
		- eValue : List of e-value float values
		- fastaFiles : List of fasta files of each monomeric PDB structures
		- pdbFiles : List of PDB Files of each monomerics structures
		- workDir : The path to the file where data and models will be generate
		- hhbin : Path to the HHSuite binarie files

	OUTPUT : 
		None

	FUNCTION :
		This function will generate for each monomeric structure involved in the same dimeric structure a
		template Object and a query Object necessary to modelize the dimeric model. For each dimeric model
		generate It call the self.serialize method and send these data to create a file where each 
		neighboring node able to generate dimeric models for a central node is refered by:
			- Chain of central nodes
			- Chain of neighboring node
			- Pir file of the model
			- The e-values of the two nodes
			- The pdb files path for each monomeric structures involved in the same dimeric structure 
			- Fasta file of the two nodes
			- WorkDir : The path to the dimeric model generated

	This function can lunch a calculation of a contact map.

	'''

	def _launchModelisation(self, edgeObj, id_struct, chains, eValue, fastaFiles, pdbFiles, workDir, hhbin):
		print "Modelisation"
		
		contactMapData = {
			"queries" : [],
			"templates" : [],
			"hhAlignFiles" : []
		}

		data = { edgeObj.v1.query : { edgeObj.v2.query : { id_struct : { } }}}

		for fasta, pdb, evalue, chain in zip(fastaFiles, pdbFiles, eValue, chains):
			
			parser = PDB.Parser()
			template = parser.load(file=pdb)
			query = qr.Query(fastaFile = fasta)
		
			# Lunch Threading
			pyproteins.services.utils.mkdir(workDir)
			thread = thr.hhAlign(query, [template], True, workDir, hhbin, 10)
			
			correctPath = workDir+'/'+edgeObj.v2.query+"_"+chain+"/"

			contactMapData["queries"].append(query)
			contactMapData["templates"].append(template)
			contactMapData["hhAlignFiles"].append(correctPath+"hhAlign.out")

			data[edgeObj.v1.query][edgeObj.v2.query][id_struct][chain] = { "pir" : correctPath+'Default.pir', "modele": workDir , "structRef" : pdb, "query" : fasta, "e-value" : str(evalue) }
		
		self._serialize(data, '/Users/mbachir/Desktop/fastaGraph/Modelisation/'+self.graphObj.queryCenter+'/'+edgeObj.v1.query+'_'+edgeObj.v2.query+'/'+id_struct+'.json')

		####### CONTACT MAP OPERATIONS ########
		#print "Check Contact Map"
		#thr.checkContactMap(contactMapData["queries"], contactMapData["templates"], contactMapData["hhAlignFiles"])
		#######################################


	'''
	
		DESERIALIZE

	INPUT : 
		- file : The JSON file contaning serialize data

	OUTPUT : 
		- A generator which raise the pdbs and queries monomeric structure paths with the pir file

	FUNCTION:
		Deserialize a JSON file

	'''

	def _deserialize(self, file):
		with open (path, 'r') as file:
			data = json.load(file)
        	for nodeC in data:
        		for nodeN in data[nodeC]:
        			for model in data[nodeC][nodeN]:
        				pdbs = [data[nodeC][nodeN][model]['pdb_monom_R'], data[nodeC][nodeN][model]['pdb_monom_L']]
        				queries = [data[nodeC][nodeN][model]['query_monom_R'], data[nodeC][nodeN][model]['query_monom_L']]
        				workDir = data[nodeC][nodeN][model]['pir']

        				yield (pdbs, queries, workDir)

    '''
			FIND BLAST
		
		INPUT : 
			- nameFile : The name of the file containing the blast output
			- path : The path to the directory containing folders with blast output inside

		OUTPUT:
			- Blast output file

		FUNCTION:
			This method will recursively search throught the tree all .blast file format and
			raise the path of each one.
    '''

	def _findBlast(self, nameFile, path):
		if os.path.isdir(path):
			for root, subdirs, filenames in os.walk(path):
				for filename in fnmatch.filter(filenames, nameFile+'.blast'):
					yield os.path.join(root, filename)
		else:
			if path.endswith('.blast'):
				yield path

	'''
			COVERAGE PROTEIN

		INPUT:
			nodeList : List of nodes 
			e-value :  float e-value

		OUTPUT : 
			- A matplotlib graph 

		FUNCTION:
			It will show a graph of cumulative models generated by the nodes for a neighbor graph
	'''

	def _coverageProtein(self, nodeList, evalue):
		
		dic = {}

		for protein in nodeList:
			for monomer in protein.monomer:
				if float(monomer[3]) <= evalue:
					if protein not in dic:
						dic[protein] = []
					dic[protein].append(monomer)
		print len(dic)
		x_axis = np.linspace(1, len(dic), len(dic))
		y_axis = [len(v) for k, v in dic.iteritems()]
		#for k, v in dic.iteritems():
		#	print k.query, len(v)


		plt.hist([y_axis], normed=1, cumulative=True, bins=150)
		plt.show()


	'''
			GRAPH MODEL PER VERTICE

		INPUT:
			edgeObj : The edge Object

		OUTPUT : 
			- A matplotlib graph 

		FUNCTION:
			It will show a graph of models generated by the neighboring nodes for a neighbor graph
	'''


	def _graphModelPerVertice(self, edgeObj):

		#Count number of model for given protein (central node)
		neighboors = []
		models = []

		for neighbor in self.modelForEdge[edgeObj.v1.query]:
			neighboors.append(neighbor)
			models.append(self.modelForEdge[edgeObj.v1.query][neighbor][0])

		x_axis = np.linspace(1, len(neighboors), len(neighboors))

		print "Graph x: Nombre de modeles y: Arretes"
		plt.figure(figsize=(15,8))
		plt.hist2d(x_axis, models, bins=60)
		plt.ylim((min(models),max(models)))

		plt.xticks(x_axis, neighboors)
		plt.show()

	'''
			GRAPH MODEL PER PROT

		INPUT:
			None

		OUTPUT : 
			- A matplotlib graph 

		FUNCTION:
			It will show a graph of models generated by a neighbor graph
	'''

	def graphModelPerProt(self):

		proteinList = []
		models = []
		count = 0

		#Count number of model for given protein (central node)
		for protein in self.modelForEdge:
			for neighbor in self.modelForEdge[protein]:
				if self.modelForEdge[protein][neighbor][0] > 0:
					count += self.modelForEdge[protein][neighbor][0]
			if count > 0:
				proteinList.append(protein)
				models.append(count)
				count = 0

		x_axis = np.linspace(1, len(proteinList), len(proteinList))
		print "Graph x: Noeuds centrals y: Nombre de modeles"

		plt.hist2d(x_axis, models, bins=20, clip_on=True)
		#plt.ylim((min(models),max(models)))

		plt.xticks(x_axis, proteinList)
		plt.show()

	'''
			BLAST SPLITTER

		INPUT:
			data : Corresponding of the data containing in a blast file 

		OUTPUT : 
			- A String data

		FUNCTION:
			Read the all the data contained in a blast file and split it whene the parser
			meets a ">" corresponding to an hit and return all the data in string format between
			the first ">" met and the next ">" 
	'''

	def _blastSplitter(self, data, separator=lambda x: x.startswith('>')):
		buff = []
		separatorOn = False
		for line in data:
			if separator(line):
				separatorOn = True
				if buff:
					yield ''.join(buff)
					buff[:] = []

			if separatorOn:
				buff.append(line)

		yield ''.join(buff)

	'''
			CLASS : EDGE

		INPUT:
			v1 : First vertice of edge
			v2 : Second vertice of edge

		ARGUMENTS:
			- hitsV1 : All hits find in blast output during a BLAST operation throught a givien DB for the first vertice
			- hitsV2 : All hits find in blast output during a BLAST operation throught a givien DB for the second vertice
	'''

class Edge(object):
	def __init__(self, v1, v2):
		self.v1 = v1
		self.v2 = v2
		self.hitsV1 = []
		self.hitsV2 = []

	'''
			CLASS : NODE

		INPUT:
			query : Uniprot ID

		ARGUMENTS:
			- monomer : All hits as monomer find in blast output during a BLAST operation throught a givien DB for node
	'''

class Node(object):
	def __init__(self, query):
		self.query = query
		self.monomer = []
