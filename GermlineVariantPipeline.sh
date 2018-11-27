#!/usr/bin/bash
## Asia Mitchell 
## UASAGE:  ./align_call.sh file.fastq(full_path) patientName sampleName sequencingPlatform(ie Illumina)
export PATH=$PATH:/opt/installed/:/usr/bin/ # Add path to bioinformatic software to PATH

fq=$1
name=$2
sample=$3
plat=$4

################
## Setup temp directory, all output files will be located here
tmp_dir='/home/users/admitch/ProjectName/'
################



echo $fq
outbam=$tmp_dir$name'-'$sample'.bam'


### REFERENCE & ANNOTATION DATA
ref='RefUCSC/g1k/human_g1k_v38.fasta'
indel='1000g_indels.vcf.gz'
thousndg='1000G_phase1.snps.high_confidence.hg38.vcf.gz'
omni1000g='1000G_omni2.5.hg38.sites.vcf.gz'
dbsnp='dbsnp_138.b38.vcf.gz'
exac='ExAC.r0.3.1.sites.vep.vcf.gz'
hapmap='hapmap_3.3.hg38.sites.vcf.gz'

### SET UP READ GROUP
readgrp='@RG\tID:'$name'\tSM:'$sample'\tPL:'$plat'\tLB:'$name

################
## STEP 1: Align FASTQ using BWA 
## FASTQ to SAM
################
if [ $? -eq 0 ]
 	then
 	
## If read length is > 75bps, use BWA mem
## If paired-end reads, give the parameter: -Map (remove 'p' for single-end reads)
		/home/exacloud/lustre1/SpellmanLab/mitcheas/User/Tools/bwa-0.7.17/bwa mem -Map -R $readgrp $ref $fq > $outsam
		
## If read length is < 75bps and reads are single-end, use BWA samse
		#/home/exacloud/lustre1/SpellmanLab/mitcheas/User/Tools/bwa-0.7.17/bwa samse -n 100 -r $readgrp $ref - $fq > $outsam
 	else
 		echo 'Exit at FASTQ to SAM:'$?
 		exit
fi
################
## Sort SAM & Convert SAM to BAM
## SAM must be sorted by chromosome/contig position.
## Convert SAM to a BAM file. BAM files are compressed and required by subsequent tools in
## the pipeline.
## Once the final BAM is created, we can delete the SAM file.
###############
if [ $? -eq 0 ]
 	then
 		out=$tmp_dir$name'.sorted'
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar SortSam \
 			-I $outsam \
 			-O $out'.bam' \
			--VALIDATION_STRINGENCY LENIENT \
 			-SO coordinate 
 	else
 		echo 'Exit at Sort SAM, SAM to BAM:'$?
 		exit
fi
################
## STEP 2: MAP AND MARK PCR DUPLICATES	
## Don't remove PCR duplicates from final BAM to preserve complete data in case of future 
## with updated pipeline.
################
if [ $? -eq 0 ]
	then
		infile=$out
		out=$outbam
		metric=$out'.metrics'
		
		rm $outsam
		
		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar MarkDuplicates \
			--ASSUME_SORTED true \
			--REMOVE_DUPLICATES false \
			--VALIDATION_STRINGENCY LENIENT \
			-I $infile \
			-O $out \
			-M $metric
	else
		echo 'Exit at MAP AND MARK DUPLICATES:'$?
		exit
fi
###############
## Build BAM index after each new BAM even if previous tool created BAM Index outfile.
## This prevents any time stamp discrepancy - sometimes the BAM index writing completes before 
## the BAM file because the index file is much smaller and quicker to write. This can cause
## time stamp discrepancies and break the pipeline.
###############
if [ $? -eq 0 ]
	then
		infile=$out
		out=$out'.bai'
		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar BuildBamIndex \
			-I $infile \
			--VALIDATION_STRINGENCY LENIENT \
			-O $out
	else
		echo 'Exit at MAP AND MARK DUPLICATES - BuildIndex:'$?
		exit
fi

#################
## CALCULATE COVERAGE
## Determine genome-wide coverage - you can also limit coverage calculation to target area
## using the '-L' parameter and a .bed file listing target regions.
#################
if [ $? -eq 0 ]
	then
		out=$infile'.depthofcov'
		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar DepthOfCoverage \
			-R $ref \
			-o $out \
			-I $infile \
			-L /home/exacloud/lustre1/users/mitcheas/VHL/LowPassWGS/chrom_arms.bed \
			--omitDepthOutputAtEachBase
	else
		echo 'Exit at Calculate coverage:' $?
		exit
fi

#################
## BASE QUALITY SCORE RECALIBRATION
## Detects systematic errors from the sequencing instrument and adjusts the quality score
## for bases with such errors. 
## Later, you may want to look at the recalibration tables across many sequence runs to
## to determine if there are any quality concerns. Use "AnalyzeCovariate" tool to perform
## this comparison. 
#################

# 1st - Create a table of loci in known-sites vcfs (loci where variation is expected) and
# of covariate values (covariates include read group, reported quality score, nucleotide 
# context, etc...).
if [ $? -eq 0 ]
	then
		out=$infile'.recal_data.table'
		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar BaseRecalibrator \
   			-I $infile \
   			-R $ref \
   			--known-sites $dbsnp \
   			--known-sites $exac \
   			-O $out
   	else
		echo 'Exit at BaseRecalibrator:' $?
		exit
fi

# 2nd - Use recalibration table to correct for any bias and print new BAM.
if [ $? -eq 0 ]
	then
		out=$infile'.recal'
		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar ApplyBQSR \
   			-I $infile \
   			-R $ref \
   			--bqsr-recal-file $out \
   			-O $out'.bam'
   	else
		echo 'Exit at ApplyBQSR:' $?
		exit
fi

# 3rd - Build BAM index.

if [ $? -eq 0 ]
	then
		infile=$out'.bam'
		out=$out'.bai'
		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar BuildBamIndex \
			-I $infile \
			--VALIDATION_STRINGENCY LENIENT \
			-O $out
	else
		echo 'Exit at BSQSR - BuildIndex:'$?
		exit
fi

#################
## GERMLINE VARIANT CALLING
## This process is memory intensive and gives the parameter "-Xm4g" to set the maximum
## heap size to 4 gigs so HaplotypeCaller has sufficient memory.
## The 'ERC' parameter is the mode for emitting reference confidence scores. The GVCF mode
## compresses the output from each base to variant blocks. This mode requires the 
## the 'variant_index_type' and 'variant_index_parameter' to be set as below.
## The 'G' parameter describes the annotation applied to the variant. 'AS_Standard' is 
## allele specific annotation which is helpful for rare variant discovery and evaluates
## alleles independently at multi-allelic sites. 
## The 'L' parameter restricts variant calling to regions listed in bed file.
##
## In addition to the sample VCF, HaplotypeCaller outputs an intermediate GVCF file which
## can be used later for efficient joint genotyping across many samples.
#################

# 1st - Create VCF for single sample and the GVCF for later cohort genotyping.
if [ $? -eq 0 ]
	then
		infile=$infile'.bam'
		out=$infile'.vcf.gz'
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
 			--java-options '-Xmx4g' HaplotypeCaller  \
   			-R $ref \
   			-I $infile'.bam' \
   			-O $out \
   			-L 'target_regions.bed' \
   			-ERC GVCF \
   			-variant_index_type LINEAR \
   			-variant_index_parameter 128000 \
   			-G Standard \
   			-G AS_Standard
	else
		echo 'Exit at HaplotypeCaller:'$?
		exit
fi

# 2nd - Import single sample VCFs to Genome database to improve efficiency of joint
#		genotyping across cohort. GVCF files for each sample are listed in the sample map
#		file 'cohort.sample_map'. 'L' must be provided and is limited to a single region.

if [ $? -eq 0 ]
	then
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
 			--java-options '-Xmx4g -Xms4g' GenomicsDBImport  \
       		--genomicsdb-workspace-path $tmp_dir \
       		--batch-size 50 \
       		-L 'chr1:1-22000000' \
      		--sample-name-map 'cohort.sample_map' \
       		--reader-threads 5
	else
		echo 'Exit at HaplotypeCaller:'$?
		exit
fi

# 3rd - Genotype across a cohort.

if [ $? -eq 0 ]
	then
		out=$tmp_dir'cohort.g.vcf.gz'
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
 			--java-options '-Xmx4g' GenotypeGVCFs  \
   			-R $ref \
   			-I gendb://$tmp_dir \
   			-O $out 
	else
		echo 'Exit at GenotypeGVCFs'$?
		exit
fi

#################
## VARIANT RECALIBRATION
## Build a model to use for recalibrating variants to apply filters. Step one builds the 
## model based on information from known polymorphic sites, provided as resources. 
## Step two applies to model and filters variants. 
#################

# 1st - Build recalibration model. 'AS' parameter performs in allele specific mode because
#		VCF was produce in AS mode. Resources given are truth and false-sets - priors set 
#		here will be altered depending on sequencing protocol (WGS, WES or target capture).
# 		Prior likelihoods are confidence that the data is reliable as a truth set.
#		'an' parameter describes the annotations to provide in the final VCF. This will 
# 		only recalibrate SNPs. Indel recalibration can be performed separately using indel
#		resource files (ie $indel).

if [ $? -eq 0 ]
	then
		out=$infile'.AS.recalVQSR'
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar VariantRecalibrator \
   			-R $ref \
   			-V $infile'.vcf.gz' \
   			-AS \
   			--resource hapmap,known=false,training=true,truth=true,prior=15.0:$hapmap \
   			--resource omni,known=false,training=true,truth=false,prior=12.0:$omni1000g \
   			--resource 1000G,known=false,training=true,truth=false,prior=10.0:$thousndg \
   			--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$dbsnp \
   			-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   			-mode SNP \
   			--recal-file $infile'.AS.recal' \
   			--tranches-file $infile'.AS.tranches' \
   			--rscript-file $infile'.AS.R' 
	else
		echo 'Exit at VariantRecalibrator'$?
		exit
fi

# 2nd - Apply recalibration model and filter variants based on model. This is only for SNPs
# 		and must be performed in 'indel' mode for indels if desired. The tranches threshold
#		can be altered to allow for more permissive filtering or to reduce false positives.

if [ $? -eq 0 ]
	then
		out=$infile'.AS.recalVQSR'
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar ApplyVQSR \
   			-R $ref \
   			-V $infile'.vcf.gz' \
   			-AS \
   			-O $infile'recal.vcf.gz' \
   			-AS \
   			--ts_filter_level 99.0 \
   			--tranches-file $infile'.AS.tranches' \
   			--recal-file $infile'.AS.recal' \
   			-mode SNP 
	else
		echo 'Exit at ApplyVQSR'$?
		exit
fi

################	
## After reviewing final VCF and metrics files - you will want to remove intermediate files
################
rm $outsam
rm $outbam
rm $outbam.bai
rm $tmp_dir$name'.sorted.bam'
rm $tmp_dir$name'.sorted.bam.bai'
rm $infile'.recal_data.table'
rm $tmp_dir$name.target_intervals.list 
rm $tmp_dir$name.metrics
rm $infile'.AS.recal'
rm $infile'.AS.tranches'
rm $infile'.AS.R' 

#################
## VARIANT ANNOTATION & EVALUATION
## Your VCF of candidate variants need to be annotated by gene/coding region, mutation
## type (nonsynonymous, missense, etct...), and functional impact.
## The tool, Funcotator provides these annotations in VCF format. 
## The VCF can then be filtered based on parameters now annotated in VCF.
################

if [ $? -eq 0 ]
	then
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar Funcotator \
   			-R $ref \
   			-V $infile'.vcf.gz' \
   			-O $infile'.ann.vcf.gz' \
   			--data-sources-path dataSourcesFolder/ \
   			--ref-version hg38
	else
		echo 'Exit at Funcotator'$?
		exit
fi

# Filter variants - Set parameters as needed depending on interest/context. The example here
# selects for variants that are heterozygous, have a variant quality score of at least 30 
# (phred), a depth of at least 75 reads, and reside in gene YFG or ABC.

if [ $? -eq 0 ]
	then
 		java -jar /opt/installed/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar VariantFiltration \
   			-R $ref \
   			-V $infile'.ann.vcf.gz' \
   			-O $infile'filt.ann.vcf.gz' \
  			--filterExpression "isHet == 1 && QUAL > 30.0 && DP == 75 && GENE_NAME == YFG || 
  				GENE_NAME == ABC" \
   			--filterName "my_filters"
	else
		echo 'Exit at VariantFiltration'$?
		exit
fi