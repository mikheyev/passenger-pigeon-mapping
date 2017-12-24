from snakemake.utils import R
SRRData = {"34.3.23.2" : ["SRR5414915"], "40360" : ["SRR5414914"], "BMNH1149": ["SRR1303448", "SRR1303453", "SRR1303457"], "BMNH794" : ["SRR1303525", "SRR1303521", "SRR1303519", "SRR1303516"],  "AMNH_DOT_14025" : ["SRR5420140"], "BTP2013" : ["SRR5420139", "SRR5420138", "SRR5420137","SRR5420136"]}

outDir = "data" # symlinked to /work

localrules: combineJellyfish

rule trim:
	input:
		read1 = outDir + "/reads/{sample}_1.fastq.gz",
		read2 = outDir + "/reads/{sample}_2.fastq.gz"
	output: 
		read1 = outDir + "/readsTrimmed/{sample}_1.fastq.gz",
		read2 = outDir + "/readsTrimmed/{sample}_2.fastq.gz",
		collapsed = outDir + "/readsTrimmed/{sample}_collapsed.fastq"
	params: trimDir = outDir + "/readsTrimmed/"
	threads: 12
	resources: time=8*60, mem=20  #time in minutes, mem in gigabytes
	shell:
		"""
		module load adapterremoval/2.2.2
		AdapterRemoval --file1 {input.read1} --file2 {input.read2} \
		--output1 {output.read1} --output2 {output.read2}  \
		--collapse --minalignmentlength 15 --outputcollapsed {output.collapsed} \
		--basename {params.trimDir}/{wildcards.sample} \
		--trimwindows 10 --minquality 20 \
		--minlength 35 --threads {threads} 
		"""

rule jellyfish:
	input: rules.trim.output.read1, rules.trim.output.read2,
	output: outDir + "/reads/jellyfish/{sample}.txt"
	threads: 12
	resources: mem=100, time=2880
	shell:
		"""
		module load jellyfish/2.2.6
		# jellyfish can swamp home directory
		cd {outDir}/reads/jellyfish
		jellyfish count -t 12 -C -m 19 -s 5G  -o 40360_out <(zcat {input})
		jellyfish histo -o 40360.txt 40360_out && rm 40360_out
		"""

rule combineJellyfish:
	input: expand(outDir + "/reads/jellyfish/{sample}.txt", sample = SRRData.keys())
	output: outDir + "/reads/jellyish/jellyfish.pdf"
	run:
		# print("""
		# require(tidyverse)
		# pdf(file = {output}, height = 2 , width = 4)
		# read_csv(pipe("for i in *txt ; do awk -v name=$(basename $i .txt) -v OFS=, \'{{print name,$1,$2}}\' $i; done "), col_names=c("id","freq","count")) %>% mutate(species=ifelse(name == "AMNH_DOT_14025" | name == "BTP2013", "bt", "pp")) %>% ggplot(aes(freq, count, color=name, shape=species))+geom_point()+scale_y_log10()+xlim(0,200)
		# dev.off()
		# """)
		R("""
		require(tidyverse)
		pdf(file = {output}, height = 2 , width = 4)
		read_csv(pipe("for i in {input} ; do awk -v name=$(basename $i .txt) -v OFS=, \'{{print name,$1,$2}}\' $i; done "), col_names=c("id","freq","count")) %>% 
		mutate(species=ifelse(name == "AMNH_DOT_14025" | name == "BTP2013", "bt", "pp")) %>% 
		ggplot(aes(freq, count, color=name, shape=species))+geom_point()+scale_y_log10()+xlim(0,200)
		dev.off()
		""")

rule all:
	input: expand("data/readsTrimmed/{sample}_collapsed.fastq", sample = SRRData.keys())
