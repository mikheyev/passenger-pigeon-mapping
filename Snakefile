from snakemake.utils import R
from collections import defaultdict
from Bio import SeqIO

SRRData = {"34.3.23.2" : ["SRR5414915"], "40360" : ["SRR5414914"], "BMNH1149": ["SRR1303448", "SRR1303453", "SRR1303457"], "BMNH794" : ["SRR1303525", "SRR1303521", "SRR1303519", "SRR1303516"],  "AMNH_DOT_14025" : ["SRR5420140"], "BTP2013" : ["SRR5420139", "SRR5420138", "SRR5420137","SRR5420136"]}
outDir = "/work/MikheyevU/sasha/passenger-pigeon-mapping/data" 
SCRATCH = "/work/scratch/sasha"
clivBWAindex = outDir + "/../ref/GCA_000337935.2_Cliv_2.1_genomic.fna"

localrules: combineJellyfish, filterDisco, discoSplit

rule trim:
	input:
		read1 = outDir + "/reads/{sample}_1.fastq.gz",
		read2 = outDir + "/reads/{sample}_2.fastq.gz"
	output: 
		read1 = outDir + "/readsTrimmed/{sample}_1.fastq.gz",
		read2 = outDir + "/readsTrimmed/{sample}_2.fastq.gz",
	params: trimDir = outDir + "/readsTrimmed"
	threads: 12
	resources: time=8*60, mem=20  #time in minutes, mem in gigabytes
	shell:
		"""
		module load adapterremoval/2.2.2
		AdapterRemoval --gzip --file1 {input.read1} --file2 {input.read2} \
		--output1 {output.read1} --output2 {output.read2}  \
		--trimwindows 10 --minquality 20 \
		--basename {params.trimDir}/{wildcards.sample} \
		--minlength 35 --threads {threads} 
		"""

rule jellyfish:
	input: rules.trim.output
	output: outDir + "/readsTrimmed/jellyfish/{sample}_}{k}.txt" # fixme
	threads: 12
	params: jellydir = outDir + "/readsTrimmed/jellyfish", kmerFile = lambda wildcards: outDir + "/readsTrimmed/jellyfish/{wildcards.sample}_{wildcards.k}_out"
	resources: mem=100, time=2880
	shell:
		"""
		module load jellyfish/2.2.6
		# jellyfish can swamp home directory
		cd {params.jellydir}
		jellyfish count -t 12 -U 100 -C -m {wildcards.k} -s 5G  -o {params.kmerFile} <(zcat {input})
		jellyfish histo -o {output} {params.kmerFile} && rm {params.kmerFile}
		"""

rule fastqc:
	input: outDir + "/readsTrimmed/{sample}_{R}.fastq.gz"
	output: outDir + "/readsTrimmed/{sample}_{R}_fastqc.html"
	threads: 12
	resources: mem=10, time=60*5
	shell:
		"""
		module load fastqc/0.11.5
		fastqc -t {threads} --noextract -f fastq {input}
		"""

rule discoSNP:
	input: expand(outDir + "/readsTrimmed/{sample}_{R}.fastq.gz", sample = SRRData.keys(), R = [1,2])
	output:
		vcf = outDir + "/discoSNP/discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf",
		fasta = outDir + "/discoSNP/discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa"
	threads: 12
	resources: mem=200, time=60*24*5
	params: discodir = outDir + "/discoSNP", samples = " ".join(SRRData.keys()), inputs = " ".join(map(lambda x: x + ".txt", SRRData.keys()))
	shell:
		""" 
		module load DiscoSnp/2.3.0
		cd {params.discodir}
		for i in {params.samples}; do
			ls -1 {outDir}/readsTrimmed/$i*fastq.gz > $i.txt
		done
		echo {params.inputs} | tr " " "\\n"  > fof.txt
		run_discoSnp++.sh -u -l -r fof.txt
		"""

# filter out sites with at most one missing value from the fasta, and output them as tab-delimited
rule filterDisco:
	input: rules.discoSNP.output
	output: outDir + "/discoSNP/filtered.txt"
	threads: 1
	run:
		keep = set([])
		with open(input[0]) as vcf:
			for line in vcf:
				if line[0] != "#":
					missing = 0
					for rec in line.split("\t")[9:]:						
						if rec[:3] == "./.":
							missing += 1
						if missing > 1:
							break
					else:
						keep.add(line.split("\t")[0].split("_")[-1]) #add variant id

		records = defaultdict(dict)
		with open(output[0], "w") as f:
			for seq in SeqIO.parse(open(input[1]),'fasta'):
				uniqueId = seq.id.split("|")[0].split("_")[-1]
				path = seq.id.split("|")[0].split("_")[1]
				# keep SNPs that have low missingness
				if seq.id[:3] == "SNP" and uniqueId in keep:
					records[uniqueId][path] = seq
			# verify that records are printed paired
			for seq in records:
				if "higher" in records[seq] and "lower" in records[seq]:
					f.write(">{}\t{}\t>{}\t{}\n".format(records[seq]["higher"].id, str(records[seq]["higher"].seq), records[seq]["lower"].id, str(records[seq]["lower"].seq)))

# # split the output files to speed up mapping
# rule discoSplit:
# 	input: rules.filterDisco.output
# 	output: expand(outDir + "/discoSNP/split/x{seq}",  seq = map(lambda x: str(x).zfill(2), range(100)))
# 	params: discoDir = outDir + "/discoSNP/split"
# 	shell:
# 		"""
# 		cd {params.discoDir}
# 		split -d -n l/100 {input}
# 		for i in x??; do sed -i 's/\\t/\\n/g' $i; done
# 		"""

# map kmers onto pigeon genome
# keep only uniquely mapped reads https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
rule discoMap:
	input: rules.filterDisco.output
	output: outDir + "/discoSNP/filtered.vcf"
	threads: 12
	resources: mem=20, time=60*24*2
	params: discoDir = outDir + "/discoSNP"
	shell: 
		"""
		module load bwa.icc/0.7.10 DiscoSnp/2.3.0 samtools
		cd {params.discoDir}
		sed 's/\\t/\\n/g' {input} > filtered.fa
		bwa samse {clivBWAindex} <(bwa aln -t {threads} {clivBWAindex} filtered.fa) filtered.fa | awk '$0~/^@/ {{print; next}} {{line1=$0; getline; if(match(line1""$0, "XA:Z") == 0 ) print line1"\\n"$0}}' > filtered.sam 
		run_VCF_creator.sh -f filtered.sam -o {output}
		"""

# rule discoMerge:
# 	input: expand(outDir + "/discoSNP/split/filtered.{seq}.vcf", seq = map(lambda x: str(x).zfill(2), range(100)))

#discoSnp renames 
rule renameDisco:

rule combineJellyfish:
	input: expand(outDir + "/readsTrimmed/jellyfish/{sample}_{k}.txt", sample = SRRData.keys(), k = [19, 25, 35, 45])
	output: outDir + "/readsTrimmed/jellyfish/jellyfish.pdf"
	params: indir = outDir + "/readsTrimmed/jellyfish"
	run:
		R("""
		require(tidyverse)
		p <- tibble(File = list.files("data/readsTrimmed/jellyfish", full.names = TRUE, pattern = "*.txt")) %>% extract(File, "id", "jellyfish/(.*).txt", remove=FALSE) %>% mutate(Data = lapply(File, read_table2, col_names=c("freq","count"))) %>% unnest(Data) %>%
		mutate(species=ifelse(id == "AMNH_DOT_14025" | id == "BTP2013", "bt", "pp")) %>%
		extract(File, "k", "jellyfish/.*_(??).txt", remove=FALSE) %>%
		ggplot(aes(freq, count, color=id, shape=species))+geom_point()+facet_grid(.~k)+scale_y_log10()+xlim(0,100)
		ggsave(filename = "{output}", p)
		""")

rule all:
	input: expand(outDir + "/readsTrimmed/{sample}_{R}_fastqc.html", sample = SRRData.keys(), R = [1,2]), outDir + "/readsTrimmed/jellyfish/jellyfish.pdf"
