### DRAFT README 

1. MAPPING READS
	Map the reads from the raw fastq files using the script map_reads2016.sh. To submit this script on our computing cluster we used map_reads2016.slurm and the input list of samples allSampList.txt 

2. GENOTYPE CALLING 
	A. gvcfPrepAndGeneration.sh generates GVCF files for each sample. It was submitted using wild2016_gatk.slurm and map_reads2016_sampleList.txt
	B. genotypeGVCFs.sh generates the raw VCF file. Before running this, you need to make sure that nothing went wrong with step 2A as the input is all files in the output directory that stores the GVCF files.  
