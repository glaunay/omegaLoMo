import urllib3
import os
import fnmatch
import re

import pyproteins.homology.dummyThreaders as thr
import pyproteinsExt.structure.coordinates as PDB
import pyproteins.homology.query as qr
import pyproteins.services.utils

class GraphToModel(object):
	def __init__(self, graphObj):
		self.graphObj = graphObj
		self.structureArray = []
		self.modelForEdge = {}
		self._modelGenerated = {}


	def fetchReferenceStructure(self):
		pass

	def _psiBlast(self, blastDir, iteraction, queryObj, dbFilePath, scriptPath, ourDir):

		#def psiblast(queryFilePath, dbFilePath, iteration, hhBinDir, scriptFile ,outDir):

	    string  = '#!/bin/bash\n'
	    string += hhBinDir+'blastpgp'
	    string += '-i ' + queryFilePath + ' -d ' + dbFilePath + ' '
	    string += '-o '+outDir+'/blast.out -m 7 -j '+iteration
    
	    with open(scriptFile+'_blast', "w") as f:
	        f.write(string)
	    rq = subprocess.call('sh '+scriptFile+'_blast' ,shell=True, stdout=subprocess.PIPE)

	def fastaFromNodes(self):
		http = urllib3.PoolManager()
		nodesID = set()
		edgeDict = self.graphObj.graphData

		for node in edgeDict:
			for neighboor in edgeDict[node]:
				nodesID.add(node)
				nodesID.add(neighboor)

		if not nodesID < 0:
			for ID in nodesID:
				
				ID = ID.query
				url = "http://www.uniprot.org/uniprot/"+ID+".fasta"
				r = http.request('GET', url, preload_content=False)

				if r:
					folder = "/Users/mbachir/Desktop/fastaGraph/"+self.graphObj.queryCenter+"/"+ID+"/"
					folder_sbtach = "fastaGraph/"+self.graphObj.queryCenter+"/"+ID+"/"
					os.makedirs(folder)

					# Generate sbatch file
					self.generateBatch("/home/mbachir/", folder_sbtach, ID, folder)
					
					with open(folder+ID+'.fasta', 'wb') as out:
					    while True:
					        data = r.read()
					        if not data:
					            break
					        out.write(data)
					r.release_conn()

	def generateBatch(self, homePath, workDir, fileName, localFolder):

		sbatch = "#!/bin/bash\n"
		sbatch +="#SBATCH --workdir="+homePath+workDir+"\n"
		sbatch +="#SBATCH --nodes 1 # Number of nodes\n"
		sbatch +="#SBATCH --ntasks=1 # Number of task per nodes\n"
		sbatch +="#SBATCH -p express-mobi\n"
		sbatch +="#SBATCH --qos express-mobi\n"
		sbatch +="#SBATCH -o "+homePath+workDir+"blastInterEvol.out # File to which STDOUT will be written\n"
		sbatch +="#SBATCH -e "+homePath+workDir+"blastInterEvol.err # File to which STDERR will be written\n"

		sbatch += "module load ncbi-blast/2.2.26\n"

		sbatch += "pwd\n"
		sbatch += "ls -lrtha\n"
		sbatch += 'echo "-->"$WORK_DIR\n'
		sbatch += 'echo "-->"$WORKDIR\n'
		sbatch += "cp "+fileName+".fasta  $WORKDIR\n"
		sbatch += "cd $WORKDIR\n"
		sbatch += "ls -lrtha\n"
		sbatch += "pwd\n"

		sbatch += "which blastpgp\n"

		sbatch += "sleep 30\n"



		sbatch += 'cmd="blastpgp -i $WORKDIR/'+fileName+'.fasta -j 1 -d '+homePath+'interEvol_db/multiFastaDB.fasta -o '+homePath+workDir+fileName+'.blast"\n'
		sbatch += "echo $cmd\n"
		sbatch += "$cmd\n"

		with open(localFolder+fileName+'.sbatch', "w") as f:
			f.write(sbatch)

	def homologuesInBlast(self, edgeObj, pathBlast = "/Users/mbachir/Desktop/fastaGraph/", pdb_files = "/Users/mbachir/Desktop/interevol/InterEvol/PDB/interEvol_structures"):

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
				edgeObj.hitsV1.append((verticeR, verticeR_idChain, verticeR_Chain, verticeR_Evalue, filesVerticeL, pdb_files))

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

	def findDimericTemplates(self, edgeObj, eValue, hhbin = '/Users/mbachir/Desktop/pipeline/hhsuite/bin/'):

		setPotentialDimericStruct = set([dimeric for dimeric in self._createSetOfDimeric(edgeObj, eValue)])
		path_table_dimeric_struct = '/Users/mbachir/Desktop/interevol/InterEvol/INTER70_REFINFO.table'
		workDir = '/Users/mbachir/Desktop/fastaGraph/Modelisation/'+self.graphObj.queryCenter+'/'+edgeObj.v1.query+'_'+edgeObj.v2.query
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

							if not edgeObj.v2.query in self._modelGenerated:
								self._modelGenerated[edgeObj.v2.query] = []

							if id_struct not in self._modelGenerated[edgeObj.v2.query]:

								self._launchModelisation([pathFiles_structA+".fasta", pathFiles_structB+".fasta"], 
									[pdb_files+"/"+id_struct+"_"+chaine_structA+".pdb", pdb_files+"/"+id_struct+"_"+chaine_structB+".pdb"], workDir, hhbin)
								
								self._modelGenerated[edgeObj.v2.query].append(id_struct)

								evalue_Dimeric[id_struct] = [(id_struct+"_"+chaine_structA, evalue_MonomerA), 
								(id_struct+"_"+chaine_structB, evalue_MonomerB)]
								
								model += 1

			self.modelForEdge[edgeObj.v1.query+"_"+edgeObj.v2.query] = [model, evalue_Dimeric]
	
	def _launchModelisation(self, fastaFiles, pdbFiles, workDir, hhbin):
		
		for fasta, pdb in zip(fastaFiles, pdbFiles):
			parser = PDB.Parser()
			template = parser.load(file=pdb)
			query = qr.Query(fastaFile = fasta)
		
			# Lunch Threading
			pyproteins.services.utils.mkdir(workDir)
			thread = thr.hhAlign(query, [template], True, workDir, hhbin, 1)



	def _findBlast(self, nameFile, path):
		if os.path.isdir(path):
			for root, subdirs, filenames in os.walk(path):
				for filename in fnmatch.filter(filenames, nameFile+'.blast'):
					yield os.path.join(root, filename)
		else:
			if path.endswith('.blast'):
				yield path

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

class Edge(object):
	def __init__(self, v1, v2):
		self.v1 = v1
		self.v2 = v2
		self.hitsV1 = []
		self.hitsV2 = []
