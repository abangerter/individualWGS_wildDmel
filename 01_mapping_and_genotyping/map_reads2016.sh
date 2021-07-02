#! /bin/bash


### load modules
	echo "loading modules"
	module load bwa/0.7.17
	module load samtools/1.9


### define variables
	echo "defining variables"
	sampID=${1}	
	threads=1
	trimDir="/home/ab5dr/trimmomatic/Trimmomatic-0.39"
	pearDir="/home/ab5dr/pear/pear-0.9.11-linux-x86_64/bin"
	picardDir="/home/ab5dr/picard"
	listDir="/scratch/ab5dr/wildDmel2016/fileList"
	inputDir="/scratch/ab5dr/wildDmel2016/fastq"
	interDir="/scratch/ab5dr/wildDmel2016/interBamNew"
	outputDir="/scratch/ab5dr/wildDmel2016/mappedBamNoDup"
	refDir="/scratch/ab5dr/wildDmel2016/misc"


### mapping the ~8 different fastq for one sample
	flowcell=($( ls ${inputDir}/*${sampID}*.fastq.gz | grep -oE 'HL[0-9A-Z]{1,}' | sort | uniq ))
	## map reads
		for cell in "${flowcell[@]}"; do
			echo "${cell}"
			lanes=($( ls ${inputDir}/${cell}*${sampID}* | grep -oE 's[1-7]{1,}' | sort | uniq ))
			
			for lane in "${lanes[@]}"; do
				echo "${lane}"
				# trim out nextera and index seq
					java -Xmx4g -jar ${trimDir}/trimmomatic-0.39.jar PE ${inputDir}/${cell}_${lane}_1_${sampID}*.fastq.gz ${inputDir}/${cell}_${lane}_2_${sampID}*.fastq.gz ${interDir}/${cell}_${lane}_1_${sampID}.P_trimm.fastq ${interDir}/${cell}_${lane}_1_${sampID}.U_trimm.fastq ${interDir}/${cell}_${lane}_2_${sampID}.P_trimm.fastq ${interDir}/${cell}_${lane}_2_${sampID}.U_trimm.fastq ILLUMINACLIP:${trimDir}/adapters/NexteraPE-PE.fa:2:30:10:8:true
				# first, merge overlapping reads
					${pearDir}/pear \
					-f ${interDir}/${cell}_${lane}_1_${sampID}.P_trimm.fastq \
					-r ${interDir}/${cell}_${lane}_2_${sampID}.P_trimm.fastq \
					-o ${interDir}/${cell}_${lane}_${sampID} \
					-j ${threads} -k
				# next, map to reference genome
					## assembled reads (i.e., those that were merged by PEAR)
					bwa mem -t ${threads} \
						-R "@RG\tID:${sampID};${cell};${lane};assembled\tSM:${sampID}\tPL:illumina\tPU:${sampID};${cell};${lane}" \
						${refDir}/combined.edit2.fa \
						${interDir}/${cell}_${lane}_${sampID}.assembled.fastq | \
					samtools view -Suh - | \
					samtools sort -@ ${threads} -o ${interDir}/${cell}_${lane}_${sampID}.assembled.sort.bam
					samtools index ${interDir}/${cell}_${lane}_${sampID}.assembled.sort.bam
					## unassembled reads
					bwa mem -t ${threads} \
						-R "@RG\tID:${sampID};${cell};${lane};unassembled\tSM:${sampID}\tPL:illumina\tPU:${sampID};${cell};${lane}" \
						${refDir}/combined.edit2.fa \
						${interDir}/${cell}_${lane}_${sampID}.unassembled.forward.fastq \
						${interDir}/${cell}_${lane}_${sampID}.unassembled.reverse.fastq | \
					samtools view -Suh - | \
					samtools sort -@ ${threads} -o ${interDir}/${cell}_${lane}_${sampID}.unassembled.sort.bam
					samtools index ${interDir}/${cell}_${lane}_${sampID}.unassembled.sort.bam
			done
		done

			
### merge bam files for a sample
	echo "merging files"
	## assembled 
		ls ${interDir}/*${sampID}.assembled.sort.bam > ${listDir}/${sampID}.assembled.txt
		samtools merge -b ${listDir}/${sampID}.assembled.txt ${interDir}/${sampID}.assembled.merge.bam
		samtools index ${interDir}/${sampID}.assembled.merge.bam
	## unassembled
		ls ${interDir}/*${sampID}.unassembled.sort.bam > ${listDir}/${sampID}.unassembled.txt
		samtools merge -b ${listDir}/${sampID}.unassembled.txt ${interDir}/${sampID}.unassembled.merge.bam
		samtools index ${interDir}/${sampID}.unassembled.merge.bam


### mark duplicates on the merged file 
	echo "marking duplicates"
	## assembled 
    	java -jar ${picardDir}/picard.jar MarkDuplicates INPUT=${interDir}/${sampID}.assembled.merge.bam OUTPUT=${outputDir}/${sampID}.assembled.merge.noDup.bam METRICS_FILE=${outputDir}/${sampID}.assembled.merge.noDup.metrics REMOVE_DUPLICATES=TRUE CREATE_INDEX=true
   ## unassembled
     	java -jar ${picardDir}/picard.jar MarkDuplicates INPUT=${interDir}/${sampID}.unassembled.merge.bam OUTPUT=${outputDir}/${sampID}.unassembled.merge.noDup.bam METRICS_FILE=${outputDir}/${sampID}.unassembled.merge.noDup.metrics REMOVE_DUPLICATES=TRUE CREATE_INDEX=true



