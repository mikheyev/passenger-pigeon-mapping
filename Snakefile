from snakemake.utils import R
from collections import defaultdict
from Bio import SeqIO

SRRData = {"34.3.23.2" : ["SRR5414915"], "40360" : ["SRR5414914"], "BMNH1149": ["SRR1303448", "SRR1303453", "SRR1303457"], "BMNH794" : ["SRR1303525", "SRR1303521", "SRR1303519", "SRR1303516"],  "AMNH_DOT_14025" : ["SRR5420140"], "BTP2013" : ["SRR5420139", "SRR5420138", "SRR5420137","SRR5420136"]}
QCgood = ["34.3.23.2", "40360", "AMNH_DOT_14025", "BTP2013"] # eliminate BMNH samples, which have poor kmer profiles
outDir = "/work/MikheyevU/sasha/passenger-pigeon-mapping/data" 
SCRATCH = "/work/scratch/sasha"
clivBWAindex = outDir + "/../ref/GCA_000337935.2_Cliv_2.1_genomic.fna" # rock dove reference
btpBWAindex = outDir + "/../ref/GCA_002029285.1_NIATT_ARIZONA_genomic.fna" # band-tailed pigeon reference

localrules: combineJellyfish, discoSplit, discoMerge, discoValidate, discoValidateBT, discoFinal, discoStats, TsTv, combineStats, windowAverage

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

rule mergeReads:
	input:
		rules.trim.output
	output: outDir + "/readsTrimmed/merged/{sample}.assembled.fastq", outDir + "/readsTrimmed/merged/{sample}.unassembled.forward.fastq", outDir + "/readsTrimmed/merged/{sample}.unassembled.reverse.fastq"
	threads: 12
	params: prefix =  lambda wildcards: outDir + "/readsTrimmed/merged/" + wildcards.sample
	resources: mem = 100, time = 60*24
	shell:
		"""
		module load pear/0.9.10
		pear -f {input[0]} -r {input[1]} -j {threads} -o {params.prefix} -n 30
		"""

rule pear:
	input: expand(outDir + "/readsTrimmed/merged/{sample}.assembled.fastq", sample = SRRData.keys())

rule jellyfish:
	input: rules.mergeReads.output
	output: outDir + "/readsTrimmed/jellyfish/{sample}_{k}.txt" # fixme
	threads: 12
	params: jellydir = outDir + "/readsTrimmed/jellyfish", kmerFile = lambda wildcards: outDir + "/readsTrimmed/jellyfish/" + wildcards.sample + "_" + wildcards.k + "_out"
	resources: mem=100, time=2880
	shell:
		"""
		module load jellyfish/2.2.6
		# jellyfish can swamp home directory
		cd {params.jellydir}
		jellyfish count -t 12 -U 100 -C -m {wildcards.k} -s 5G  -o {params.kmerFile} <(cat {input})
		jellyfish histo -o {output} {params.kmerFile} && rm {params.kmerFile}
		"""

rule combineJellyfish:
	input: expand(outDir + "/readsTrimmed/jellyfish/{sample}_{k}.txt", sample = SRRData.keys(), k = [19, 25, 35, 45])
	output: outDir + "/readsTrimmed/jellyfish/jellyfish.pdf"
	params: indir = outDir + "/readsTrimmed/jellyfish"
	run:
		R("""
		require(tidyverse)
		p <- tibble(File = list.files("data/readsTrimmed/jellyfish", full.names = TRUE, pattern = "*.txt")) %>% extract(File, "id", "jellyfish/(.*)_\\\\d+.txt", remove=FALSE) %>% mutate(Data = lapply(File, read_table2, col_names=c("freq","count"))) %>% unnest(Data) %>%
		mutate(species=ifelse(id == "AMNH_DOT_14025" | id == "BTP2013", "bt", "pp")) %>% extract(File, "k", "jellyfish/.*_(\\\\d+).txt", remove=FALSE) %>%
		ggplot(aes(freq, count, color=id, shape=species))+geom_point()+facet_grid(k~., scales = "free_y")+scale_y_log10()+xlim(0,50)+theme_bw()+theme(legend.position="bottom")
		ggsave(filename = "{output}", p)
		""")

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

#to do use merged reads QCgood
rule discoSNP:
	input: expand(outDir + "/readsTrimmed/merged/{sample}.{R}.fastq", sample = QCgood, R = ["assembled", "unassembled.forward", "unassembled.reverse"])
	output:
		vcf = outDir + "/discoSNP/discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf",
		fasta = outDir + "/discoSNP/discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa"
	threads: 20
	resources: mem=500, time=60*24*2
	params: discodir = outDir + "/discoSNP", samples = " ".join(QCgood), inputs = " ".join(map(lambda x: x + ".txt", QCgood)), sampleDir = outDir + "/readsTrimmed/merged"
	shell:
		""" 
		module load DiscoSnp/2.3.0
		cd {params.discodir}
		for i in {params.samples}; do
			ls -1 {params.sampleDir}/$i*assembled.fastq -1 {params.sampleDir}/$i*unassembled*fastq > $i.txt
		done
		echo {params.inputs} | tr " " "\\n"  > fof.txt
		run_discoSnp++.sh -u -t -l -r fof.txt
		"""

# filter out sites with missing values from the fasta, and output them as tab-delimited
rule filterDisco:
	input: rules.discoSNP.output
	output: outDir + "/discoSNP/filtered.txt"
	threads: 1
	resources: mem=80, time=60*24
	run:
		keep = set([])
		with open(input[0]) as vcf:
			for line in vcf:
				if line[0] != "#":
					count = 0
					for rec in line.split("\t")[9:]:						
						if rec[:3] == "./.":
							count += 1
						if count > 0:   #legacy
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


# map kmers onto Rock pigeon genome
# keep only uniquely mapped reads https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
rule discoMap:
	input: rules.filterDisco.output
	output: outDir + "/discoSNP/disco.vcf"
	threads: 12
	resources: mem=20, time=60*24*2
	params: discoDir = outDir + "/discoSNP"
	shell: 
		"""
		module load miniconda DiscoSnp/2.3.0 samtools
		cd {params.discoDir}
		sed 's/\\t/\\n/g' {input} > filtered.fa
		bwa samse {clivBWAindex} <(bwa aln -t {threads} {clivBWAindex} filtered.fa) filtered.fa |  awk '$0~/^@/ {{print; next}} {{line1=$0; getline; if(match(line1""$0, "XA:Z") == 0 ) print line1"\\n"$0}}' > filtered.sam 
		run_VCF_creator.sh -f filtered.sam -o {output}
		"""


rule discoMapBT:
	input: rules.discoSNP.output[1]
	output: outDir + "/discoSNP/discoBT.vcf"
	threads: 12
	resources: mem=20, time=60*24*2
	params: discoDir = outDir + "/discoSNP"
	shell: 
		"""
		module load miniconda DiscoSnp/2.3.0 samtools
		cd {params.discoDir}
		bwa samse {btpBWAindex} <(bwa aln -t {threads} {btpBWAindex} {input}) {input} |  awk '$0~/^@/ {{print; next}} {{line1=$0; getline; if(match(line1""$0, "XA:Z") == 0 ) print line1"\\n"$0}}' | grep -v "INDEL" > filteredBT.sam 
		run_VCF_creator.sh -f filteredBT.sam -o {output}
		"""

def validate(refFile, inFileName, outFileName ):
	# validates results of discoSNP mapping and creates a filtered vcf
	seqs = SeqIO.to_dict(SeqIO.parse(refFile, "fasta"))
	outfile = open(outFileName, "w")
	badSites = set()
	sites = {}
	with open(inFileName) as vcf:
		for line in vcf:
			if line[0] == "#":
				outfile.write(line)
				continue
			line = line.split("\t")
			if "_path_" in line[0]:
				#remove unmapped variants
				continue
			siteId = (line[0], line[1])
			if siteId in sites:
				# if there is a duplicate, filter site
				badSites.add(siteId)
				continue
			else:
				sites[siteId] = line
			if (int(line[1]) - 1) > len(seqs[line[0]]) or str(seqs[line[0]][int(line[1])-1]) != line[3]:
				# reference doesn't match
				badSites.add(siteId)
	for site in sites:
		if site not in badSites:
			outfile.write("\t".join(sites[site]))

# verify that reference positions match, and eliminate any sites that have multiple
rule discoValidate:
	input: rules.discoMap.output
	output: outDir + "/discoSNP/disco.validated.vcf"
	run:
		validate(clivBWAindex, input[0], output[0])

rule discoValidateBT:
	input: rules.discoMapBT.output
	output: outDir + "/discoSNP/discoBT.validated.vcf"
	run:
		validate(btpBWAindex, input[0], output[0])

#filter results to remove sited more than 3*sd above the mean coverage, and anything with other than a PASS flag
rule discoFinal:
	input: rules.discoValidate.output
	output: outDir + "/discoSNP/disco.final.vcf"
	shell:
		"""
		module load vcftools
		function mean_sd () {{
				awk '{{x[NR]=$0; s+=$0; n++}} END{{a=s/n; for (i in x){{ss += (x[i]-a)^2}} sd = sqrt(ss/n); print s/n,sd}}'
				}}
		stats=( $(awk '$1!~/^#/ {{sum=0; for(i=10;i<=NF;i++) {{split($i,a,":"); sum+=a[2]}}; print sum/4}}' {input} | mean_sd ) )
		echo DP mean: ${{stats[0]}} DP stdev: ${{stats[1]}}
		cutoff=$(python -c "from math import ceil; print(int(ceil(${{stats[0]}} +  3 * ${{stats[1]}})))")
		vcftools --vcf {input} --remove-filtered-all --max-meanDP $cutoff --recode --stdout > {output} && [[ -s {output} ]]
		"""

#filter results to remove sited more than 3*sd above the mean coverage, and anything with other than a PASS flag
rule discoFinalBT:
	input: rules.discoValidateBT.output
	output: outDir + "/discoSNP/discoBT.final.vcf"
	shell:
		"""
		module load vcftools
		function mean_sd () {{
				awk '{{x[NR]=$0; s+=$0; n++}} END{{a=s/n; for (i in x){{ss += (x[i]-a)^2}} sd = sqrt(ss/n); print s/n,sd}}'
				}}
		stats=( $(cat {input} | vcftools --vcf - --max-missing 1 --recode --stdout | awk '$1!~/^#/ {{sum=0; for(i=10;i<=NF;i++) {{split($i,a,":"); sum+=a[2]}}; print sum/4}}' | mean_sd ) )
		echo DP mean: ${{stats[0]}} DP stdev: ${{stats[1]}}
		cutoff=$(python -c "from math import ceil; print(int(ceil(${{stats[0]}} +  3 * ${{stats[1]}})))")
		vcftools --vcf {input} --remove-filtered-all --max-missing 1 --max-meanDP $cutoff --recode --stdout > {output} && [[ -s {output} ]]
		"""

#snpEff using a pre-built database from gff3
rule snpEff:
	input: rules.discoFinal.output
	output: outDir + "/discoSNP/disco.annotated.vcf"
	resources: mem=10, time = 60*24
	shell:
		"""
		java -Xmx7g -jar /apps/free/snpeff/4.3g/snpEff.jar -c /apps/unit/MikheyevU/sasha/snpEff4/snpEff.config -no-utr -no-upstream -no-intron -no-intergenic -no-downstream cliv_2.1 {input} > {output} && [[ -s {output} ]]
		"""

rule snpEffBT:
	input: rules.discoFinalBT.output
	output: outDir + "/discoSNP/discoBT.annotated.vcf"
	resources: mem=10, time = 60*24
	shell:
		"""
		java -Xmx7g -jar /apps/free/snpeff/4.3g/snpEff.jar -c /apps/unit/MikheyevU/sasha/snpEff4/snpEff.config -no-utr -no-upstream -no-intron -no-intergenic -no-downstream btp {input} > {output} && [[ -s {output} ]]
		"""


# compute various stats from discoSNP data
#G1,G2 are PP 34.3.23.2 and 40360
#G3,G4 are BT AMNH_DOT_14025 and BTP2013
rule discoStats:
	input: rules.snpEff.output
	output: 
		PPMK = outDir + "/reports/PP_MK.txt",
		BPMK = outDir + "/reports/BP_MK.txt",
		PPpi = outDir + "/reports/PP_pi.txt",
		BPpi = outDir + "/reports/BP_pi.txt",
		PPfreq = outDir + "/reports/PP_freq.txt",
		BPfreq = outDir + "/reports/BP_freq.txt"
	shell:
		"""
		module load vcftools
		vcftools --vcf {input} --recode-INFO-all --indv G1 --indv G2 --non-ref-af .01 --recode --stdout | \
		awk -v OFS="\\t" '$8~/synonymous_variant/ {{print $1,$2,"S"}} $8~/missense_variant/ {{print $1,$2,"N"}}' > {output.PPMK} && [[ -s {output.PPMK} ]]

		vcftools --vcf {input} --recode-INFO-all --indv G3 --indv G4 --non-ref-af .01 --recode --stdout | \
		awk -v OFS="\\t" '$8~/synonymous_variant/ {{print $1,$2,"S"}} $8~/missense_variant/ {{print $1,$2,"N"}}' > {output.BPMK} && [[ -s {output.BPMK} ]]

		vcftools --vcf {input} --indv G1 --indv G2 --site-pi --stdout > {output.PPpi} && [[ -s {output.PPpi} ]]

		vcftools --vcf {input} --indv G3 --indv G4 --site-pi --stdout > {output.BPpi} && [[ -s {output.BPpi} ]]

		vcftools --vcf {input} --indv G1 --indv G2 --freq2 --stdout | cut -f-5 > {output.PPfreq} && [[ -s {output.PPfreq} ]]

		vcftools --vcf {input} --indv G3 --indv G4 --freq2 --stdout | cut -f-5 > {output.BPfreq} && [[ -s {output.BPfreq} ]]
		"""

# combine the various selection stats
rule combineStats:
	input:
		PPMK = outDir + "/reports/PP_MK.txt",
		BPMK = outDir + "/reports/BP_MK.txt",
		PPpi = outDir + "/reports/PP_pi.txt",
		BPpi = outDir + "/reports/BP_pi.txt",
		PPfreq = outDir + "/reports/PP_freq.txt",
		BPfreq = outDir + "/reports/BP_freq.txt"
	output:
		pi = outDir + "/reports/pi.csv",
		PPMK = outDir + "/reports/PP_MK.csv",
		BPMK = outDir + "/reports/BP_MK.csv"
	run:
		R("""
			library(tidyverse)
			left_join(read_tsv("{input.PPpi}"), read_tsv("{input.BPpi}"), by = c("CHROM" = "CHROM", "POS" = "POS")) %>% rename(PPpi = PI.x , BPpi = PI.y) %>% write_csv( path="{output.pi}")
			left_join(read_tsv("{input.PPMK}", col_names=c("CHROM", "POS", "type")),
                  read_tsv("{input.PPfreq}", 
                           col_names=c("CHROM", "POS", "N_ALLELES", "N_CHR","freq"), skip = 1),
                   by =  c("CHROM" = "CHROM", "POS" = "POS")) %>% select(-(4:5)) %>%
    			mutate(poly = ifelse(freq > 0 & freq < 1, "P", "D")) %>% 
    			mutate(variant = paste0(poly, type)) %>% write_csv(path = "{output.PPMK}")
    		left_join(read_tsv("{input.BPMK}", col_names=c("CHROM", "POS", "type")),
                  read_tsv("{input.BPfreq}", 
                           col_names=c("CHROM", "POS", "N_ALLELES", "N_CHR","freq"), skip = 1),
                   by =  c("CHROM" = "CHROM", "POS" = "POS")) %>% select(-(4:5)) %>%
			    mutate(poly = ifelse(freq > 0 & freq < 1, "P", "D")) %>% 
    			mutate(variant = paste0(poly, type)) %>% write_csv(path = "{output.BPMK}")
    		""")

# compare TsTv for reads mapped to rock pigeon vs band-tailed pigeon
rule TsTv:
	input: 
		RP = rules.discoFinal.output,
		BT = rules.discoValidateBT.output
	output: 
		RP = outDir + "/reports/TsTv_RP.txt",
		BT = outDir + "/reports/TsTv_BT.txt"
	shell:
		"""
		module load vcftools
		# collect Ts and Tv data for SNPs that have some difference from ancestor
		paste \
		<(vcftools --vcf {input.RP} --non-ref-af .01 --TsTv-summary --indv G1 --indv G2  --stdout) \
		<(vcftools --vcf {input.RP} --non-ref-af .01 --TsTv-summary --indv G3 --indv G4  --stdout | cut -f2 ) | \
		sed '1cMODEL\\tPP\\tBT' > {output.RP} && [[ -s {output.RP} ]]

		# calculate TsTv ratios for SNPs that all mapped to the BT reference, without missing data
		paste \
		<(vcftools --vcf {input.BT} --max-missing 1 --recode --stdout | vcftools --vcf - --non-ref-af .01 --TsTv-summary --indv G1 --indv G2 --stdout) \
		<(vcftools --vcf {input.BT} --max-missing 1 --recode --stdout | vcftools --vcf -  --non-ref-af .01 --TsTv-summary --indv G3 --indv G4 --stdout | cut -f2) | \
		sed '1cMODEL\\tPP\\tBT' > {output.BT} && [[ -s {output.BT} ]]
		"""

rule all:
	input: outDir + "/reports/TsTv_RP.txt", outDir + "/reports/pi.csv"
