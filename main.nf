#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Help message
def helpMessage = """
Usage: nextflow run main.nf [options]

Options:
  --input_dir       Input directory containing the input bam files (required)
  --output_dir      Output directory (required)
  --ref             Path to the Hg38 reference genome (required)
  --annotationsDir  Path to AnnotSV annotation directory (optional), otherwise annotation of merged SV calls will be skipped
  --HPO_terms       HPO terms for SV annotation using SvAnna from individual SV callers [cuteSV, sniffles2, svim] (optional)
                    :HPO terms should be provided as txt file as follows "e.g. --phenotype-term HP:0001249 --phenotype-term HP:0001250", otherwise SvAnna annotation will be skipped
  --STR_regions     STR regions for quantification using NanoRepeat (optional)
                    :Regions to genotype STRs should be provided as a TAB separated file "e.g. chrX\t148500638\t148500683\tCCG", otherwise STR quantification will be skipped
  --str             Flag to trigger STR genotyping using nanoRepeat (requires --STR_regions)
  --consensus       Flag to trigger all processes except nanoRepeat
  
Other flags:
  --help           Print this help message

"""

// Parse command-line arguments
params.help = false
params.str = false
params.consensus = false

if (params.help || !params.input_dir || !params.output_dir || !params.ref) {
    println helpMessage
    exit 0
}

// Check if output directory exists
def outputDir = file(params.output_dir)
if (outputDir.exists()) {
    println "Output directory exists: ${params.output_dir}"
} else {
    println "Output directory does not exist. Creating one: ${params.output_dir}"
    outputDir.mkdirs()
}

// worfkflow params
params.input_dir = ""
params.output_dir = ""
params.ref = ""
params.HPO_terms = ""
params.STR_regions = ""
params.annotationsDir = ""
params.tr_bed_path = "human_GRCh38_no_alt_analysis_set.trf.bed"
tr_bed = file("${workflow.projectDir}/${params.tr_bed_path}")

// To run SvAnna annotation, HPO terms should be provided with --HPO_terms as follows "e.g. --phenotype-term HP:0001249 --phenotype-term HP:0001250"; otherwise no SvAnna annotation will be skipped
// To run NanoRepeat STR quantification, a file containing locus (tab seperated) to genotype STRs should be provided with --STR_regions as follows "e.g. chrX\t148500638\t148500683\tCCG

input_ch = Channel.fromPath(params.input_dir + '/*.bam', checkIfExists: true)
    .map { [it.baseName, file(it), file(it + '.bai')] }

    
/* ------------------------- cuteSV ------------------------- */

// To install cuteSV either do: {pip install cuteSV} or {conda install -c bioconda cutesv} or {git clone https://github.com/tjiangHIT/cuteSV.git && cd cuteSV/ && python setup.py install}  

process cutesv {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(name), file(bam), file(bai)

    output:
    path("${bam.baseName}_cutesv.filtered.vcf"), emit: cutesv_vcf

    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    cuteSV ${bam} ${params.ref} ${bam.baseName}_cutesv.vcf ./ \
        -S ${bam.baseName} \
        -s 4 \
        -l 50 \
        -L 100000 \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --genotype \
        --gt_round 500 \
        --threads ${task.cpus} &&
    
    bcftools sort ${bam.baseName}_cutesv.vcf -o ${bam.baseName}_cutesv.sorted.vcf
    bcftools view -i 'FILTER="PASS"' ${bam.baseName}_cutesv.sorted.vcf -o ${bam.baseName}_cutesv.filtered.vcf
    """
}

/* ------------------------- SvAnna ------------------------- */

//SvAnna - Efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. For installation see: {https://github.com/TheJacksonLaboratory/SvAnna}
//SvAnna annotation doesn't work well with the combiSV output. Seems to ignore all deletions. Better to run SvAnna on each SV caller output independently

process annotate_cutesv_vcf {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/SvAnna", mode: 'copy'
    cpus 4
    memory '8 GB'
    
    input:
    path(cutesv_vcf)
    path(HPO_terms)
    
    output:
    path("${cutesv_vcf.baseName}.SvAnna.*")

    script:
    """
    cp ${params.HPO_terms} phenotype_terms.txt

    source /opt/conda/bin/activate lrs_analysis
    java -jar /opt/svanna-cli-1.0.4/svanna-cli-1.0.4.jar prioritize \
        --data-directory=/opt/svanna-data/ \
        --output-format html,tsv,vcf \
        --out-dir ./ \
        --prefix ${cutesv_vcf.baseName}.SvAnna \
        --n-threads ${task.cpus} \
        --uncompressed-output \
        --report-top-variants 50 \
        --vcf ${cutesv_vcf} `cat phenotype_terms.txt`
    """
}

/* ------------------------- SVIM ------------------------- */

// To install svim either do: {pip install svim} or {conda install --channel bioconda svim} or {git clone https://github.com/eldariont/svim.git && cd svim && pip install .}  

process svim {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls", mode: 'copy'
    cpus 8
    memory '16 GB'
    
    input:
    tuple val(name), file(bam), file(bai)

    output:
    path("${bam.baseName}_svim.filtered.vcf"), emit: svim_vcf
    
    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    svim alignment ./ ${bam} ${params.ref} \
        --sample ${bam.baseName} \
        --min_sv_size 50 \
        --max_sv_size 100000 \
        --types DEL,INS,INV,DUP:TANDEM,DUP:INT \
        --minimum_depth 4  &&
    
    bcftools sort variants.vcf -o sorted.variants.vcf
    bcftools view -i 'QUAL >= 5' sorted.variants.vcf -o sorted.QUAL.variants.vcf
    bcftools view -i 'FILTER="PASS"' sorted.QUAL.variants.vcf -o ${bam.baseName}_svim.filtered.vcf
    """
}

/* ------------------------- SvAnna ------------------------- */

//SvAnna - Efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. For installation see: {https://github.com/TheJacksonLaboratory/SvAnna}
//SvAnna annotation doesn't work well with the combiSV output. Seems to ignore all deletions. Better to run SvAnna on each sv caller output independently

process annotate_svim_vcf {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/SvAnna", mode: 'copy'
    cpus 4
    memory '8 GB'
    
    input:
    path(svim_vcf)
    path(HPO_terms)
    
    output:
    path("${svim_vcf.baseName}.SvAnna.*")

    script:
    """
    cp ${params.HPO_terms} phenotype_terms.txt

    source /opt/conda/bin/activate lrs_analysis
    java -jar /opt/svanna-cli-1.0.4/svanna-cli-1.0.4.jar prioritize \
        --data-directory=/opt/svanna-data/ \
        --output-format html,tsv,vcf \
        --out-dir ./ \
        --prefix ${svim_vcf.baseName}.SvAnna \
        --n-threads ${task.cpus} \
        --uncompressed-output \
        --report-top-variants 50 \
        --vcf ${svim_vcf} `cat phenotype_terms.txt`
    """
}

/* ------------------------- Sniffles2 ------------------------- */

// To install svim either do: {conda install sniffles=2.0}

process sniffles2 {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls", mode: 'copy'
    cpus 8
    memory '16 GB'
    
    input:
    tuple val(name), file(bam), file(bai)

    output:
    path("${bam.baseName}_sniffles.filtered.vcf"), emit: sniffles_vcf
    
    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    sniffles \
        --threads ${task.cpus} \
        --sample-id ${bam.baseName} \
        --reference ${params.ref} \
        --minsvlen 50 \
        --symbolic \
        --cluster-merge-pos 150 \
        --tandem-repeats ${tr_bed} \
        --input ${bam} \
        --vcf ${bam.baseName}_sniffles.vcf &&
        
    sed '/.:0:0:0:NULL/d' ${bam.baseName}_sniffles.vcf > tmp.vcf
    mv tmp.vcf ${bam.baseName}_sniffles.vcf
    bcftools view -i 'FILTER="PASS"' ${bam.baseName}_sniffles.vcf -o ${bam.baseName}_sniffles.filtered.vcf
    """
}

/* ------------------------- SvAnna ------------------------- */

//SvAnna - Efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. For installation see: {https://github.com/TheJacksonLaboratory/SvAnna}
//SvAnna annotation doesn't work well with the combiSV output. Seems to ignore all deletions. Better to run SvAnna on each SV caller output independently

process annotate_sniffles2_vcf {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/SvAnna", mode: 'copy'
    cpus 4
    memory '8 GB'
    
    input:
    path(sniffles_vcf)
    path(HPO_terms)

    output:
    path("${sniffles_vcf.baseName}.SvAnna.*")
    
    script:
    """
    cp ${params.HPO_terms} phenotype_terms.txt

    source /opt/conda/bin/activate lrs_analysis
    java -jar /opt/svanna-cli-1.0.4/svanna-cli-1.0.4.jar prioritize \
        --data-directory=/opt/svanna-data/ \
        --output-format html,tsv,vcf \
        --out-dir ./ \
        --prefix ${sniffles_vcf.baseName}.SvAnna \
        --n-threads ${task.cpus} \
        --uncompressed-output \
        --report-top-variants 50 \
        --vcf ${sniffles_vcf} `cat phenotype_terms.txt`
    """
}

/* ------------------------- combiSV ------------------------- */

//combiSV combines sample SV calls from sniffles, svim, and cuteSV to generate a high-confident call SV VCF per sample. {see: https://github.com/ndierckx/combiSV.git}

process mergeCalls {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/CombiSV", mode: 'copy'
    cpus 2
    memory '4 GB'
    
    input:
    tuple val(name), file(cutesv_vcf), file(svim_vcf), file(sniffles_vcf)
    
    output:
    path("${name}_combiSV.merged.consensus.vcf"), emit: consensus_vcf

    script:
    """
    perl /opt/combiSV/combiSV2.3.pl \
        -sniffles ${sniffles_vcf} \
        -cutesv ${cutesv_vcf} \
        -svim ${svim_vcf} \
        -o ${name}_combiSV.merged.consensus.vcf
    
    """
}

/* ------------------------- AnnotSV ------------------------- */

//AnnotSV - An integrated tool for Structural Variations annotation and ranking. For installation see: {https://github.com/lgmgeo/AnnotSV}
//Annotates the output of consensus calls from combiSV and emits tsv,logs from Annotsv

process annotate_mergedCalls {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/AnnotSV", mode: 'copy'
    cpus 12
    memory '16 GB'
    
    input:
    path(merged_vcf)
    
    output:
    path("${merged_vcf.baseName}.Annotsv.tsv")

    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    AnnotSV \
        -annotationsDir ${params.annotationsDir} \
        -SVinputFile ${merged_vcf} \
        -outputFile ./${merged_vcf.baseName}.Annotsv.tsv \
        -snvIndelPASS 1 \
        -SVminSize 50 \
        -genomeBuild GRCh38 \
        > ${merged_vcf.baseName}.Annotsv.log
    """
}

/* ------------------------- NanoRepeat ------------------------- */

//NanoRepeat: quantification of Short Tandem Repeats (STRs) from long-read sequencing data
//To install and set up NanoRepeat, see https://github.com/WGLab/NanoRepeat#installation 

process nanoRepeat {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/STRs/", mode: 'copy'
    cpus 4
    memory '12 GB'

    input:
    tuple val(name), file(bam), file(bai)

    output:
    path("${bam.baseName}.details")
    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    nanoRepeat.py \
        --input ${bam} \
        --type bam \
        --ref_fasta ${params.ref} \
        --repeat_region_bed ${params.STR_regions} \
        --num_cpu ${task.cpus} \
        -d ont_q20 \
        -o ${bam.baseName}
    """
}

workflow {
    if (params.str && !params.consensus) {
        // Run STR analysis only
        if (params.STR_regions != "" && params.STR_regions != null) {
            STR_analysis_results = nanoRepeat(input_ch)
        } else {
            println("--STR_regions not provided: Skipping quantification of STRs using NanoRepeat")
        }
    } else {
        // Run SV calling processes
        cutesv_calling_results = cutesv(input_ch)
        svim_calling_results = svim(input_ch)
        sniffles_calling_results = sniffles2(input_ch)

        // Run SvAnna annotation if HPO terms are provided
        if (params.HPO_terms != "" && params.HPO_terms != null) {
            cutesv_vcf_annotation_results = annotate_cutesv_vcf(cutesv_calling_results, params.HPO_terms)
            svim_vcf_annotation_results = annotate_svim_vcf(svim_calling_results, params.HPO_terms)
            sniffles_vcf_annotation_results = annotate_sniffles2_vcf(sniffles_calling_results, params.HPO_terms)
        } else {
            println("--HPO_terms not provided: Skipping SvAnna annotation of VCFs from individual SV callers [CuteSV, Sniffles2, SVIM]")
        }

        cutesv_calls = cutesv_calling_results.cutesv_vcf.map { file -> [file.baseName.split('_')[0], file] }
        svim_calls = svim_calling_results.svim_vcf.map { file -> [file.baseName.split('_')[0], file] }
        sniffles_calls = sniffles_calling_results.sniffles_vcf.map { file -> [file.baseName.split('_')[0], file] }

        // Collect the output and merge the VCF files for each sample separately to run consensus SV merging in combiSV process
        merged_vcfs = cutesv_calls.join(svim_calls, by: 0).join(sniffles_calls, by: 0)
        combisv_results = mergeCalls(merged_vcfs)

        // Pass the output from the combiSV process as input to the Annotsv_annotation process
        if (params.annotationsDir != "" && params.annotationsDir != null) {
            Annotsv_annotation_results = annotate_mergedCalls(combisv_results.consensus_vcf)
        } else {
            println("--annotationsDir not provided: Skipping AnnotSV annotation of consensus SVs")
        }

        // Run STR analysis if --str is also provided
        if (params.str && params.STR_regions != "" && params.STR_regions != null) {
            STR_analysis_results = nanoRepeat(input_ch)
        } else if (params.str) {
            println("--STR_regions not provided: Skipping quantification of STRs using NanoRepeat")
        }
    }
}