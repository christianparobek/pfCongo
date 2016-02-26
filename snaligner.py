## A Snakefile
## Adapted for aligning the AMA1-antibody adaped
## Pf3D7 strains from Dr. Sheetij Dutta's group
## Begun 10 March 2015


##########################################################################################


######## Some Definitions ##########
workdir: '/proj/meshnick/users/ChristianP/congoHRP2/'
readWD = '/proj/meshnick/users/ChristianP/congoHRP2/'

REF = '/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta'
GATK = '/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar'
PICARD = '/nas02/apps/picard-1.88/picard-tools-1.88'
TMPDIR = '/netscr/prchrist/tmp_for_picard/'

SAMPLES, = glob_wildcards('/proj/meshnick/users/ChristianP/congoHRP2/symlinks/{sample}_R1.fastq.gz')


##########################################################################################


####### Target #######
rule all:
#	input: expand('aln/{sample}.bam', sample = SAMPLES)
#	input: expand('aln/{sample}.dedup.bam', sample = SAMPLES)
#	input: expand('coverage/{sample}.cov10', sample = SAMPLES)
#	input: 'names/bamnames.list'
#	input: 'coverage/coverage.txt'
	input: 'coverage/cov_plot.pdf'
#	input: 'names/bamnames.list'
#	input: 'variants/05xAT100%.intervals'
#	input: 'coverage/cov_plot.pdf'




rule remove_non_entries:
	input: 'variants/split_indivs_HC/{sample}_HC.vcf'
	output: 'variants/split_indivs_HC/{sample}_HC.sans0.vcf'
	shell: 'grep -vP "PL\t0" {input} | grep -vP "\tGT\t." > {output}'

rule split_vcf:
	input: 'variants/all10_HC.pass.vcf'
	output: 'variants/split_indivs_HC/{sample}_HC.vcf'
	shell: 'java -jar {GATK} -T SelectVariants \
		-R {REF} --variant {input} \
		-sn {wildcards.sample} \
		-o {output}'

rule select_records:
	input: 'variants/all10_HC.qual.vcf'
	output: 'variants/all10_HC.pass.vcf'
	shell: 'java -jar {GATK} -T SelectVariants \
	-R {REF} -V {input} -o {output}\
	-select "vc.isNotFiltered()" \
	-restrictAllelesTo BIALLELIC'

rule variant_filtration:
	input: vcf = 'variants/all10_HC.vcf', intervals = 'variants/all10_HC_05xAT100%.intervals'
	output: 'variants/all10_HC.qual.vcf'
	shell: 'java -jar {GATK} -T VariantFiltration \
	-R {REF} -V {input.vcf} -o {output} \
	-L {input.intervals} \
	--filterExpression "QD < 10.0" --filterName "QD" \
	--filterExpression "MQ < 50.0" --filterName "MQ" \
	--filterExpression "FS > 60.0" --filterName "FS" \
	--filterExpression "MQRankSum < -5.0" --filterName "MQRankSum" \
	--filterExpression "ReadPosRankSum < -5.0" --filterName "ReadPosRankSum" \
	--logging_level ERROR'

rule filter_by_depth:
	input: 'variants/all10_HC.vcf'
	output: 'variants/all10_HC_05xAT100%.intervals'
	shell: 'java -jar {GATK} -T CoveredByNSamplesSites \
		-R {REF} -minCov 05 -percentage 0.99999 \
		-V {input} -out {output}'

rule unified_genotyper_together :
	input: bams = 'names/bamnames.list'
	output: 'variants/all10_UG.vcf'
	threads: 1
	shell: 'java -jar {GATK} -T UnifiedGenotyper \
		-R {REF} -I {input.bams} \
		-nt {threads} -ploidy 1 -o {output}'
		# all_chrs.intervals includes only chrs and mito

rule make_bamnames_list:
	input: 'names/sample.list'
	output: 'names/bamnames.list'
	shell: 'sed s#^#aln/# names/sample.list > names/bamnames.list'

rule make_sample_list:
	input: expand('aln/{sample}.realn.bam', sample = SAMPLES)
	output: 'names/sample.list'
	shell: 'ls aln | grep .realn.bam > names/sample.list'

rule plot_coverage:
	input: 'coverage/covPlotter.r'
	output: 'coverage/cov_plot.pdf'
	shell: 'Rscript {input}'

rule make_R_cov_script:
	input: 'coverage/coverage.txt'
	output: 'coverage/covPlotter.r'
	shell: """echo '## Read in data
		coverage <- read.table("coverage/coverage.txt", header=FALSE)

		## Order data acendingly
		coverage <- coverage[with(coverage, order(V2)), ]

		## Add a third column with numbers in it
		coverage$V5 <- 1:length(coverage$V1)

		##Plot the coverage graph
		pdf(file="coverage/cov_plot.pdf", width = 4, height = 4)
		plot(coverage$V2 ~ coverage$V5, axes=FALSE, xlab="", ylab="Frac Genome Covered", ylim=c(0,1), col="black", pch=20)
		points(coverage$V3 ~ coverage$V5, col="grey", pch=20)
		points(coverage$V4 ~ coverage$V5, col="red", pch=20)
		axis(1, at=1:length(coverage$V1), labels=coverage$V1, cex.axis=.7, las=3, cex = 0.5)
		axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0))
		legend(2, 0.5, c("at 5x", "at 10x", "at 25x"), pch=20, col=c("black","grey","red"))
		dev.off()' > coverage/covPlotter.r
		"""

rule digest_coverage:
	input: expand('coverage/data/{sample}.{cov}', sample = SAMPLES, cov = 'cov05 cov10 cov25'.split())
	output: 'coverage/coverage.txt'
	shell: 'for name in `ls coverage/data/ | grep cov05 | sed "s/\.cov..//"`; \
		do \
		cov05=$(tail -1 coverage/data/$name.cov05 | cut -f 5); \
		cov10=$(tail -1 coverage/data/$name.cov10 | cut -f 5); \
		cov25=$(tail -1 coverage/data/$name.cov25 | cut -f 5); \
		echo -e $name"\t"$cov05"\t"$cov10"\t"$cov25 >> coverage/coverage.txt; \
		done'

rule calculate_25x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov25'
	shell: 'bedtools genomecov \
		-ibam {input} -max 25 | grep genome \
		> {output}'

rule calculate_10x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov10'
	shell: 'bedtools genomecov \
		-ibam {input} -max 10 | grep genome \
		> {output}'

rule calculate_05x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov05'
	shell: 'bedtools genomecov \
		-ibam {input} -max 5 | grep genome \
		> {output}'

rule realn_indels:
	input: bam = 'aln/{sample}.dedup.bam', targets = 'aln/{sample}.realigner.intervals', 
	output: 'aln/{sample}.realn.bam'
	shell: 'java -jar {GATK} -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-targetIntervals {input.targets} \
		-o {output}' 
		# all_chrs.intervals includes just chrs and mito
		# -fixMisencodedQuals must be added for SRA data

rule find_indels:
	input: bam = 'aln/{sample}.dedup.bam', index = 'aln/{sample}.dedup.bai'
	output: 'aln/{sample}.realigner.intervals'
	shell: 'java -jar {GATK} -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-o {output}'
		# all_chrs.intervals includes just  chrs and mito

rule index_dedups: 
	input: 'aln/{sample}.dedup.bam'
	output: 'aln/{sample}.dedup.bai'
	shell: 'java -jar {PICARD}/BuildBamIndex.jar INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'

rule mark_dups:
	input: 'aln/{sample}.sorted.bam'
	output:'aln/{sample}.dedup.bam','aln/{sample}.dedup.metrics'
	shell: 'java -jar {PICARD}/MarkDuplicates.jar \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=TRUE \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'

rule sort_bam:
	input: 'aln/{sample}.bam'
	output: 'aln/{sample}.sorted.bam'
	shell: 'samtools sort {input} -o {output}'

rule fastq_to_bam:
	input: 'symlinks/{sample}_R1.fastq.gz', 'symlinks/{sample}_R2.fastq.gz'
	output: 'aln/{sample}.bam'
	shell: 'bwa mem {REF} {readWD}{input[0]} {readWD}{input[1]} \
		-R "@RG\tID:bwa\tPL:illumina\tLB:{wildcards.sample}\tSM:{wildcards.sample[0]}{wildcards.sample[1]}{wildcards.sample[2]}{wildcards.sample[3]}{wildcards.sample[4]}" \
		-M -t 8 -v 2 -A 2 -L 15 -U 9 -T 75 \
		-k 19 -w 100 -d 100 -r 1.5 -c 10000 \
		-B 4 -O 6 -E 1 | samtools view -Sb - > {output}'
		# calling the @RG ID: 'bwa' because this resolves a clash with @PG ID
