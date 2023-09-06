# Michelle Webb
# December 19, 2022
# Extract just germline sample SNP sites for further analysis

## the samples listed in this step are the original HC VCF:
for i in `ls *.vcf`; do j=`echo $i | cut -d- -f1`; echo ${j} ; bcftools view -s ${j} ${i} > ${j}_subset.vcf  ; done

## Remove multiallelic variants and filter for 0/1
for i in `ls justGermline/*.vcf`; do j=`echo ${i} | cut -f2 -d/ | cut -f1 -d- | cut -f1,2 -d_`; bcftools view --types snps -m 2 -M 2 -g het ${i} > ${j}_onlyHeterozygous_onlyBiAllelic.vcf ; done

## Remove alternative contigs
for i in `ls justGermline_onlyHeterozygous_onlyBiAllelic/*.vcf`; do j=`echo ${i} | cut -f2 -d/ | cut -f1 -d- | cut -f1,2 -d_`; echo ${j}; grep -v '^KI' ${i} | grep -v '^GL' > ${j}_noALTContigs.vcf; done

## Annotate - using a slurm script, use the following commands
module load usc openjdk
java -jar /scratch2/michelgw/reference/tools/snpEff/snpEff.jar -v GRCh38.105 -canon -nodownload ${vcf} > /scratch2/michelgw/reference/haplotypeCallerVCF/annotated_noQUALfilter/${sample}_annotated.vcf

## Renamed samples to match FFSlide...etc naming

## Removed alt contigs from header from all VCF, then removed chr from FFD1 vcf
cd /scratch2/michelgw/reference/haplotypeCallerVCF/annotated_noQUALfilter/
for i in `ls *vcf`; do j=`echo $i | cut -f1 -d_`; grep -v '=KI' ${i} | grep -v '=GL' > ${j}_noALTContigs_inHeader.vcf ; done
sed 's/chr//g' FFD1.vcf

## Removed M from header
for i in `ls *vcf`; do j=`echo $i | cut -f1 -d_`; grep -v '=M' ${i} > ${j}_noQUALfilter_noALTContig_noMheader.vcf ; done

## Run allele counts, filter for functional variants and QUAL afterwards
