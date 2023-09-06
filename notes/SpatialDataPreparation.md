## Spatial Data Preparation
Author: Michelle Webb         
Date: September 5, 2023       

The following example command was run for spaceranger pipeline analysis of individual samples processed with the 10X Visium platform.


**Spaceranger Pipeline**

```
/tools/spaceranger-1.1.0/spaceranger count --id=[sample] --transcriptome /scratch/working/resources/references/Homo_sapiens/10x_GRCh38-2020A/refdata-gex-GRCh38-2020-A/ --image [sample]_image.tif --slide V19N18-063 --area A1 --reorient-images --fastqs [fastq directory]
```
