#Spatial Loss of Heterozygosity Identification

The diagram below describes our bioinformatic approach for loss of heterozygosity (LOH) identification in spatial clusters.

## Contents
- Bioinformatic Approach
- Example Analysis
- Dependencies
- Contact

## Bioinformatic Approach

![image](analysisPipeline.png)

#### Exome Sequencing Germline/Tumor Analysis
Heterozygous SNP sites are first obtained through an exome bioinformatic analysis. FASTQs for germline and tumor are individually aligned to the human GRCh38 reference genome using bwa. GATK base recalibration and duplicate marking tools are applied to the resultant BAMs. GATK HaplotypeCaller is run with a list of both the germline and tumor BAMs as input. The snpEff tool annotates the joint VCF with information from dbSNP 148 (hg38). This VCF is filtered for heterozygous SNP sites. The file ExomeDataPreparation.md contains the primary commands used in this process.

#### Spatial Transcriptomics Pre-processing
Separately, the 10X Genomics spaceranger pipeline (v1.1.0) is run using spatial FASTQs as input. An example command used to start this pipeline is shown in SpatialDataPreparation.md. The two main files that proceed to the following steps of our method are the overall spatial BAM and a cluster .csv file. By default, spaceranger outputs k-means and graph cluster files. This example describes use of the default graph cluster .csv file for a sample. However, clusters can also be determined through use of alternative external tools such as Seurat or SCANPY. These clusters can be input if the file contains two comma separated columns of Barcode and Cluster.

Samples are split by cluster assignments using a python script which was originally based on a 10X Genomics utility which is now subset-bam: https://github.com/10XGenomics/subset-bam. Next, samples are indexed and secondary alignments are removed using samtools v1.9.

A custom script is run to calculate allele counts at known heterozygous SNP sites with the exome sequencing VCF file and primary alignment cluster BAMs as input. The results are filtered in R by strict criteria to obtain high quality variants. The filtered allele count .csv files are merged and converted to VCF format. This VCF is the primary input for our R package tLOH.

#### Spatial Transcriptomics LOH Identification
A sample VCF is imported in R as a dataframe. Bayes factors are calculated at each SNP site. A Hidden Markov model segmentation approach is applied in a per-chromosome, per-cluster system. Cumulative metrics across segments are evaluated to make state determinations of heterozygous, LOH, or undefined. Two plotting functions are provided for output visualization.

The main analysis functions for this process are available through Bioconductor at https://www.bioconductor.org/packages/release/bioc/html/tLOH.html. However, the most recent functions will not be available until the upcoming fall release in October. Therefore, an R script with all updated functions is available in this repository as functions.R.

## Example Analysis

An example analysis is described in the exampleLOH_analysis.html vignette. To run an example analysis from pre-procesing, the following files will be used as input:          

1. spatial transcriptomics BAM # available at NCBI GEO
2. graphclusters.csv # available at NCBI GEO
3. filtered and annotated exome VCF file # example in references directory

## Dependencies

Exome Analysis:       
bwa v17        
samtools v1.9       
gatk v4.0.10.1           
snpEff v4.3         

Spatial Pre-processing Analysis:         
spaceranger v1.1.0                        
Python 3.6.15 (all package dependencies shown in spatialProcessingEnvironment.yml)         

Spatial LOH Analysis: 
R (>= 4.2.0)      
**R packages**          
scales v1.2.1       
tidyverse v2.0.0         
ggplot2 v3.4.2         
data.table         
purrr v1.0.1        
dplyr v1.1.2         
VariantAnnotation  v.1.44.1      
GenomicRanges v.1.50.2         
naniar v1.0.0       
depmixS4 v1.5.0        
stringr v1.5.0       
stats v4.2.1        
bestNormalize v1.9.0        


## Conda environment setup

Specific python (v3.6) packages are required to run the preProcessing_tLOH.sh script. The spatial ProcessingEnvironment.yml file lists these functions. To generate a conda environment which contains these packages, the following commands can be run.

```
conda env create --name spatialEnv --file=spatialProcessingEnvironment.yml
```
Alternatively, a new conda environment can be created and each package/version listed in the .yml can be installed. For example:

```
conda create -n myEnv python=3.6
conda activate myEnv
conda config --add channels r
conda config --add channels bioconda
conda install -c bioconda
conda install pysam=0.15.4
conda install [package=version] # add package and version name
```

## Contact
**Michelle G. Webb**      
michelgw@usc.edu

## Acknowledgments
**10X Visium Spatial Gene Expression** https://www.10xgenomics.com/products/spatial-gene-expression              
**R:** R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.          
**depmixS4** Visser I, Speekenbrink M (2010). “depmixS4: An R Package for Hidden Markov Models.” Journal of Statistical Software, 36(7), 1–21. https://www.jstatsoft.org/v36/i07/.          
**scales:** Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization. R package version 1.1.1. https://CRAN.R-project.org/package=scales                          
**tidyverse:** Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686          
**ggplot2:** H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York 2016.                
**data.table:** Matt Dowle and Arun Srinivasan (2020). data.table: Extension of \`data.frame\`. R package version 1.13.0. https://CRAN.R-project.org/package=data.table                   
**purrr:** Lionel Henry and Hadley Wickham (2020). purrr: Functional Programming Tools. R package version 0.3.4. https://CRAN.R-project.org/package=purrr                           
**dplyr:** Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2020). dplyr: A Grammar of Data Manipulation. R package version 1.0.0. https://CRAN.R-project.org/package=dplyr                                 
**VariantAnnotation:** Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M (2014).
“VariantAnnotation: a Bioconductor package for exploration and annotation of
genetic variants.” _Bioinformatics_, *30*(14), 2076-2078. doi:
10.1093/bioinformatics/btu168 (URL:
https://doi.org/10.1093/bioinformatics/btu168).                 
**GenomicRanges:** Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software
  for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118.
  doi:10.1371/journal.pcbi.1003118                        
**Bayes factors** Jeffreys, Harold (1998) [1961]. The Theory of Probability(3rd ed.). 
Oxford, England. p. 432. ISBN 9780191589676.                 
**Best Normalize** Peterson RA (2021). “Finding Optimal Normalizing Transformations via bestNormalize.” The R Journal, 13(1), 310–329. doi:10.32614/RJ-2021-041.
Peterson RA, Cavanaugh JE (2020). “Ordered quantile normalization: a semiparametric transformation built for the cross-validation era.” Journal of Applied Statistics, 47(13-15), 2312-2327. doi:10.1080/02664763.2019.1630372.                  
**Naniar** https://cran.r-project.org/web/packages/naniar/index.html




