# Passenger pigeon mapping

#SRRData = {"34.3.23.2" : ["SRR5414915"], "40360" : ["SRR5414914"], "BMNH1149": ["SRR1303448", "SRR1303453", "SRR1303457"], "BMNH794" : ["SRR1303525", "SRR1303521", "SRR1303519", "SRR1303516"], 
SRRData = {"AMNH DOT 14025" : ["SRR5420140"], "BTP2013" : ["SRR5420139", "SRR5420138", "SRR5420137","SRR5420136"]}

#hed raw data for the rock dove (SRA+ SRR516969), w
flatten = lambda l: [item for sublist in l for item in sublist]
SRR = flatten(list(SRRData.values()))

outDir = "data"   # symlink to /work
refDir = "ref" 	  # symlink to /work
SCRATCH  = "/work/scratch/sasha"

localrules: all, renameDownloads

rule downloadReads:
	output: outDir + "/reads/{sample}_1.fastq.gz", outDir + "/reads/{sample}_2.fastq.gz"
	resources: mem = 3, time = 7*24*60
	shell:
		"""
		module load sratoolkit/2.8.0
		fastq-dump --skip-technical --gzip --split-files -O {outDir}/reads {wildcards.sample}
		"""

rule renameDownloads:
	input: expand(outDir + "/reads/{sample}_{R}.fastq.gz", sample = SRR, R = [1, 2])
	run:
		for sample in SRR:
			for read in SRR[sample]:
				print(read)

rule all:
	input: expand(outDir + "/reads/{sample}_{R}.fastq.gz", sample = SRR, R = [1, 2])