# Repeat expansion (STRs) and structural variant (SVs) analysis on long-read (ONT) data
`ont-research_str-and-sv-analysis` is a pipeline to call tandem repeat expansions (STRs) and structural variants (SVs) with annotation on long-read sequencing (LRS). This pipeline is inspired by my work on 
Epilepsy Research, and the fact that LRS has the potential to comprehensively identify all medically relevant genome variations, including complex SVs and STRs associated with diseases 
that are commonly missed by short-read sequencing (SRS) approaches.

The pipeline is built using [Nextflow](https://www.nextflow.io/), a bioinformatics workflow manager that enables the development of portable and reproducible workflows. 
There is a Docker container that includes all the tools/softwares required by the pipeline, which makes the results highly reproducible. 

# Repeat expansions:
Repeat expansions analysis is performed by [NanoRepeat](https://github.com/WGLab/NanoRepeat). This analysis is optional and is only triggered when regions for STR quantification are provided with the flag `--STR_regions`.
Regions to genotype STRs should be provided as a TAB separated list "e.g. chrX 148500638 148500683 CCG", otherwise STR quantification will be skipped. 

# Structural variants:
Structural variants are called by 3 pipelines; [sniffles2](https://github.com/fritzsedlazeck/Sniffles/releases), [cuteSV](https://github.com/tjiangHIT/cuteSV), and [SVIM](https://github.com/eldariont/svim), 
which call SVs based on different metrics. [combiSV](https://github.com/ndierckx/combiSV) is then applied to combine SV calls from the 3 pipelines into a superior call set per sample. This superior 
call then set for annotation by [AnnotSV](https://github.com/lgmgeo/AnnotSV). 

# Structural variant annotation:
