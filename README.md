# Repeat expansion (STRs) and structural variant (SVs) analysis on long-read (ONT) data
`ont-research_str-and-sv-analysis` is a pipeline to call tandem repeat expansions (STRs) and structural variants (SVs) with annotation on long-read sequencing (LRS) bam files. This pipeline is inspired by my work on 
Epilepsy Research, and the fact that LRS has the potential to comprehensively identify all medically relevant genome variations, including complex SVs and STRs associated with diseases 
that are commonly missed by short-read sequencing (SRS) approaches.

The pipeline is built using [Nextflow](https://www.nextflow.io/), a bioinformatics workflow manager that enables the development of portable and reproducible workflows. 
There is a Docker container that includes all the tools/softwares required by the pipeline, which makes the results highly reproducible. 

# Repeat expansions:
Repeat expansions analysis is performed by [NanoRepeat](https://github.com/WGLab/NanoRepeat). This analysis is optional and is only triggered when regions for STR quantification are provided with the flag `--STR_regions`.
Regions to genotype STRs should be provided as a TAB separated list "e.g. chrX 148500638 148500683 CCG", otherwise STR quantification will be skipped. 

# Structural variants:
Structural variants are called by 3 pipelines; [sniffles2](https://github.com/fritzsedlazeck/Sniffles/releases), [cuteSV](https://github.com/tjiangHIT/cuteSV), and [SVIM](https://github.com/eldariont/svim), 
which call SVs based on different metrics. [combiSV](https://github.com/ndierckx/combiSV) is then applied to combine SV calls from the 3 pipelines into a superior call set per sample. The resulting superior 
call set is then annotated by [AnnotSV](https://github.com/lgmgeo/AnnotSV). 

# Structural variant annotation:
Within the pipeline, structural variants can also be annotated using two annotation tools, [AnnotSV](https://github.com/lgmgeo/AnnotSV) and [SvAnna](https://github.com/TheJacksonLaboratory/SvAnna/tree/master).
SvAnna is an annotation tool that was recently developed to annotate SVs from long-read sequencing data. However, it does not annotate `DELs` from combiSV, thus it is better to annotate SVs from each caller individually using SvAnna. AnnotSV can annotate the superior call set from combiSV, thus enabling accurate annotation of clinically relevant SVs from these two sources. 

To use [SvAnna](https://github.com/TheJacksonLaboratory/SvAnna/tree/master), you need to provide HPO terms (optional). HPO terms should be provided as follows "e.g. `--phenotype-term HP:0001249 --phenotype-term HP:0001250"`, otherwise SvAnna annotation will be skipped. 

To use [AnnotSV](https://github.com/lgmgeo/AnnotSV) (required), you need to first run the `INSTALL_annotations.sh` in a directory, which downloads and installs annotation sources for AnnotSV, then provide this directory 
using this flag `--annotationsDir`. 

# Installation and Usage:
`$ git clone https://github.com/NyagaM/ont-research_str-and-sv-analysis.git` 

To view options:

`$ cd ont-research_str-and-sv-analysis`

`$ nextflow run main.nf --help`

To run the workflow, use the following command:

```bash
nextflow run main.nf -profile standard \
    --input_dir path/to/folder/with/bams \
    --output_dir path/to/results \
    --ref path/to/ref.fasta(hg38) \
    --HPO_terms path/to/HPO_terms.txt \
    --STR_regions path/to/locusIDs_for_STR_genotyping.txt \
    --annotationsDir path/to/annotsv/annotation/folder/
