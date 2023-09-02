## Exome Sequencing Data Preparation
Author: Michelle Webb        
Date: June 17, 2023 

The following commands were used in the bioinformatic processing of the exome sequencing data. Both tumor and germline samples were included in the BAMLIST variable at the HaplotypeCaller step. The output of the following step was a paired tumor and germline annotated VCF for use in the spatial transcriptomics analysis.

```
FASTQ1='example.r1.fq' 
FASTQ2='example.r2.fq' 
KNOWN='dbSNP' # version 147 
REF='Gencode_v29_GRCh38_primary_assembly.fa'
FAI='Gencode_v29_GRCh38_primary_assembly.fa.fai' 
```

**Alignment** 

```
bwa mem -R ${RGTAG} -M -t16 ${REF} ${FASTQ1} ${FASTQ2} > ${BAMPRE}.sam
samtools view -h -b -t ${FAI} ${BAMPRE}.sam -o ${BAMPRE}.bam 
samtools sort -T ${BAMPRE} ${BAMPRE}.bam -o ${BAMPRE}.bam 
samtools index ${BAMPRE}.bam 
```

**Base Recalibration** 

```
gatk --java-options "-Xmx44G" BaseRecalibrator \ 
    --reference ${REF} \ 
    --input ${BAMFILE} \ 
    --known-sites ${KNOWN} \ 
    --output ${BAMFILE}.recal_data.grp > ${BAMFILE}.recalibrateOut 


gatk --java-options "-Xmx44G" ApplyBQSR \ 
    --reference ${REF} \ 
    --input ${BAMFILE} \ 
    --output ${RECALBAM} \ 
    --bqsr-recal-file ${BAMFILE}.recal_data.grp 
```

**Mark Duplicates** 

```
gatk --java-options "-Xmx20G" MarkDuplicates \ 
    --TMP_DIR ${tempDir} \ 
    --ASSUME_SORTED=true \ 
    --VALIDATION_STRINGENCY=SILENT \ 
    --REMOVE_DUPLICATES=false \ 
    --INPUT=${BAMFILE} \ 
    --OUTPUT=${OUTPUTBAM} \ 
    --METRICS_FILE=${BAMFILE}.picStats.MarkDupMetrics \ 
    --MAX_RECORDS_IN_RAM=18000000 \ 
    --CREATE_INDEX=true 
```

**HaplotypeCaller** 

```
gatk --java-options "-Xmx44g" HaplotypeCaller \ 
    --input ${BAMLIST} \ 
    --reference ${REF} \ 
    --intervals ${CHRLIST}/Step${STEP}.list \ 
    --dbsnp ${KNOWN} \ 
    --min-base-quality-score 10 \ 
    --output ${TRK}_Step${STEP}.HC.vcf > ${TRK}_Step${STEP}.hcOut 
```

**snpEff** 

```
java -Xmx21g -jar ${SNPEFFPATH}/snpEff.jar eff \ 
    -verbose \ 
    -i vcf \ 
    -o vcf \ 
    -noLog \ 
    -stats ${summaryOut} \ 
    -noDownload \ 
    -config ${SNPEFFPATH}/snpEff.config \ 
    ${DBVERSION} \ 
    ${VCF} > $snpEffTxt 
    
java -Xmx21g -jar ${SNPEFFPATH}/snpEff.jar eff \ 
    -verbose \ 
    -i vcf \ 
    -o vcf \ 
    -noLog \ 
    -noDownload \ 
    -stats ${summaryOut} \ 
    -config ${SNPEFFPATH}/snpEff.config \ 
    ${DBVERSION} \ 
    ${VCF} > $snpEffInt 
    
java -Xmx21g -jar ${SNPEFFPATH}/SnpSift.jar annotate \ 
    ${DBSNP} \ 
    $snpEffInt > $snpEffOut 
```
