
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


#############################################



	# Input :
	# - blastDir : Blast Directory
	# - iteration : Number of iterations
	# - queryObj : Query Object made by the class Query in pyprotein
	# - hhBinDir : The path of HHSuite Binaries
	# - dbFilePath : The path to the database
	# - scriptPath : The path of the generate bash script
	# - ourDir : The path of the psiBlast Output

	# Output:
	# - Bash script to execute psi-blast

	# Function:
	# The script will generate a bash script to execute a psi_blast job throught a given database with a given
	# number of iteration in XML format 

def _psiBlast(blastDir, iteration, queryObj, hhBinDir, dbFilePath, scriptPath, ourDir):

    string  = '#!/bin/bash\n'
    string += hhBinDir+'blastpgp'
    string += '-i ' + queryFilePath + ' -d ' + dbFilePath + ' '
    string += '-o '+outDir+'/blast.out -m 7 -j '+iteration

    with open(scriptFile+'_blast', "w") as f:
        f.write(string)
    rq = subprocess.call('sh '+scriptFile+'_blast' ,shell=True, stdout=subprocess.PIPE)

	# FASTA FROM NODES
	
	# Input : No
	# Output : Folders containing FASTA files
	
	# Function :
	# Take the graph Object from the self.graphObj argument of the class
	# For each protein node of the graph It will search the corresponding
	# Fasta format sequence by It uniprot sequence inside a folder with the same name 
	# It calls in the same time the self.generateBatch() class method

def fastaFromNodes():
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

	'''
	
	FASTA FROM FILE

	Input :
	- path : Take a path file containing list a uniprot ID

	Output : Folders containing FASTA files

	Function : 
	Parse a file line by line to recover uniprot ID.
	Each Uniprot ID is used to download the corresponding Fasta File
	which is stock in a folder using the same name.
	It calls in the same time the self.generateBatch() class method
	'''

def fastaFromFile(path):
	http = urllib3.PoolManager()
	
	with open(path, "r") as file:
		for line in file:
			ID = line.replace("\n", "").split("_")[1]
			url = "http://www.uniprot.org/uniprot/"+ID+".fasta"
			r = http.request('GET', url, preload_content=False)

			if r:
				folder = "/Users/mbachir/Desktop/fastaGraph/"+ID+"/"
				folder_sbtach = "fastaGraph/"+ID+"/"
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

	'''
	GENERATE BATCH

	Input : 
	- homePath : Path of cluster home e.g : /home/mbachir/
	- workDir : Folder containing fasta files inside their correponding folders 
	- fileName : The uniprot Id used to name fasta file for a given ID
	- localFolder : Folder in the local machin containing the fasta file for a given uniprot ID

	Output :
	- A sbatch file to execute a BLAST operation on cluster
	
	Function : 
	Generate a sbtach file for a each uniprot ID that used to execute BLAST on cluster
	'''

def generateBatch(homePath, workDir, fileName, localFolder):

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