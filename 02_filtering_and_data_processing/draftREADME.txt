### DRAFT README 


### FILTER GENOTYPE FILE
The steps to filter the VCF generated with genotypeGVCFs.sh:
	A. vqsr_filtering.sh performs variant quality score recalibration, does basic filtering based on indels and repeat regions, and does some preparation for the next filtering steps.
	B. rd01_filtering.R preps the first steps of read depth based filtering
	C. rd02_filtering.sh filters the VCF file based on the calculations done in rd01_filtering.R
	D. rd03_filtering.R preps the next steps of read depth based filtering
	E. rd04_filtering.sh filters the VCF file based on the calcuations done inf rd03_filtering.R and also removes non-biallelic sites. It also submits a script that alters genotypes to missing data based on read depth filtering with scripts submitRscript.slurm that submits modifyToMissing.R, which requires files generated in rd01_filtering.R and rd03_filtering.R
	F. cnv_manta.sh detects structural variants and runs as an array using a list made in rd01_filtering.R and bam files generated during the mapping steps
	G. cnv_processing.sh takes the output from step F to make a table that contains information about some de novo structural variants (insertions, deletions, duplications)
	H. cnv_makeBedFile.R takes the table made in step G and gets list of SNPs that do not overlap with CNVs 
	I. cnv_applyFiltering.sh applies the table of CNVs made in step H and extracts them from the VCF generating the VCF file used in downstream analyses.


### GENERATE FILES FROM WHICH OTHER INFORMATION IS EXTRACTED
The script getSeqArrayGDSfile.R generates the gds file used by R package SeqArray for many analyses in this repository. 

The script getSNPRelateGDSfile.R generates the gds file used by R package SNPRelate for some analyses in this repository.


### EXTRACT BASIC INFORMATION
The script getHeterozygosityInfo.R extracts information about heterozygosity levels in the samples. This is used in last filtering steps to generate the final dataset. It will require the file generated in script getSeqArrayGDSfile.R 

The script getReadDepthInfo.R generates a table of read depth information. It will require a file generated in script getSeqArrayGDSfile.R

The script getLDprunedSNPset.R generates the list of LD-pruned SNPs. The file generated and used in our analyses can be found in the inputFiles directory. It will require a file generated by the script getSNPRelateGDSfile.R


### ANALYSES AND FIGURE MAKING
The script finalFiltering_fig01.R lays out the different final filtering steps of the dataset and generates Figure 01. It will require files generated by scripts getHeterozygosityInfo.R, getSeqArrayGDSfile.R and getReadDepthInfo.R

The script inversionTyping_fig02.R contains the method for identifying the inverted, standard, or heterozygous state of large cosmopolitan inversions in the samples and generates Figure 02. It will require files generated by scripts getSNPRelateGDSfile.R, getSeqArrayGDSfile.R, getLDprunedSNPset.R, and getHeterozygosityInfo.R. It will require the input data file v6liftOver_inv_fixed_coord.txt