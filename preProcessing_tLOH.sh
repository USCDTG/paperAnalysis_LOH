#!/bin/bash

time=`date +%d-%m-%Y-%H-%M`
echo "Start:$time"
echo ' '

conda activate spatialEnv

sample='Oligodendroglioma11'
testDir='testDirectory/'
scripts=`echo ${testDir}/scripts`
references=`echo ${testDir}/references`

echo '# Splitting BAM to cluster mini-BAMs'
python ${scripts}/splitBAMs.py --bam ${references}/${sample}_possorted_genome_bam.bam --clusters ${references}/${sample}_graphclusters.csv --output ${testDir}/outputBAM

outputBAM=`echo ${testDir}/outputBAM`

echo '# Creating output directories'
mkdir ${outputBAM}/primaryAlignments
mkdir ${testDir}/alleleCounts
mkdir ${testDir}/filteredAlleleCounts

echo '# Removing secondary alignments, indexing BAMs, and calclulating allele counts for each cluster'
echo ' '
for i in `ls ${outputBAM}/*bam`; do \

	cluster=`basename ${i}`

	INDEX=`echo ${cluster} | sed 's/[^0-9]//g'`
	echo '# Cluster' ${INDEX}
	echo ' '
 	samtools index ${i}
	samtools view -q30 -b ${i} -F 256 > ${outputBAM}/primaryAlignments/${sample}_cluster${INDEX}_primaryAlignments.bam 
	samtools index ${outputBAM}/primaryAlignments/${sample}_cluster${INDEX}_primaryAlignments.bam
	python ${scripts}/pythonAnalysis.py ${references}/${sample}_noQUALfilter_noALTContig_noM.vcf ${outputBAM}/primaryAlignments/${sample}_cluster${INDEX}_primaryAlignments.bam ${sample}_${INDEX} ${testDir}/alleleCounts
done

echo ' '
echo 'Filtering SNP sites'
Rscript ${scripts}/filter.R ${testDir}/alleleCounts
cp -v ${testDir}/alleleCounts/${sample}*filtered*.csv ${testDir}/filteredAlleleCounts
echo '# Merging allele count csv files and converting to VCF format'
echo ' '
python ${scripts}/convert_to_VCF.py ${sample} ${testDir}/filteredAlleleCounts


time=`date +%d-%m-%Y-%H-%M`
echo "End: $time"
