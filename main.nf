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

// --- MODIFIED INPUT CHANNEL ---
// Extracts sample name before the first '.' and passes it along with bam and bai files
input_ch = Channel.fromPath(params.input_dir + '/*.bam', checkIfExists: true)
    .map { bamFile ->
        // Extract sample name: part of the filename before the first dot
        def sampleName = bamFile.name.split('\\.')[0]
        // Construct the expected BAI file path
        def baiFile = file(bamFile.toString() + '.bai')
        // Check if BAI file exists (optional but recommended)
        if (!baiFile.exists()) {
            log.warn "Warning: BAI file not found for ${bamFile.name}. Expected at: ${baiFile}"
            // Decide how to handle missing BAI: error out or proceed?
            // For now, let's proceed but the process might fail later.
            // To error out: exit 1, "BAI file not found for ${bamFile.name}"
        }
        // Emit the tuple: [sampleName, bamFile, baiFile]
        [ sampleName, bamFile, baiFile ]
    }


/* ------------------------- cuteSV ------------------------- */

// To install cuteSV either do: {pip install cuteSV} or {conda install -c bioconda cutesv} or {git clone https://github.com/tjiangHIT/cuteSV.git && cd cuteSV/ && python setup.py install}

process cutesv {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/cutesv/${name}", mode: 'copy', pattern: "${name}.cutesv.filtered.vcf" // Organize output by sample name
    cpus 8
    memory '16 GB'

    input:
    tuple val(name), file(bam), file(bai) 

    output:
    tuple val(name), path("${name}.cutesv.filtered.vcf"), emit: cutesv_vcf

    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    cuteSV ${bam} ${params.ref} ${name}.cutesv.vcf ./ \
        -S ${name} \
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

    bcftools sort ${name}.cutesv.vcf -o ${name}.cutesv.sorted.vcf
    bcftools view -i 'FILTER="PASS"' ${name}.cutesv.sorted.vcf -o ${name}.cutesv.filtered.vcf
    """
}

/* ------------------------- SvAnna ------------------------- */

//SvAnna - Efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. For installation see: {https://github.com/TheJacksonLaboratory/SvAnna}
//SvAnna annotation doesn't work well with the combiSV output. Seems to ignore all deletions. Better to run SvAnna on each SV caller output independently

process annotate_cutesv_vcf {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/SvAnna/cutesv/${name}", mode: 'copy' // Organize output by sample name
    cpus 4
    memory '8 GB'

    input:
    tuple val(name), path(vcf_file)
    path(HPO_terms)

    output:
    path("${name}.cutesv.SvAnna.*")

    script:
    """
    cp ${params.HPO_terms} phenotype_terms.txt

    source /opt/conda/bin/activate lrs_analysis
    java -jar /opt/svanna-cli-1.0.4/svanna-cli-1.0.4.jar prioritize \
        --data-directory=/opt/svanna-data/ \
        --output-format html,tsv,vcf \
        --out-dir ./ \
        --prefix ${name}.cutesv.SvAnna \
        --n-threads ${task.cpus} \
        --uncompressed-output \
        --report-top-variants 50 \
        --vcf ${vcf_file} `cat phenotype_terms.txt`
    """
}

/* ------------------------- SVIM ------------------------- */

// To install svim either do: {pip install svim} or {conda install --channel bioconda svim} or {git clone https://github.com/eldariont/svim.git && cd svim && pip install .}

process svim {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/svim/${name}", mode: 'copy', pattern: "${name}.svim.filtered.vcf" 
    cpus 8
    memory '16 GB'

    input:
    tuple val(name), file(bam), file(bai) 

    output:
    tuple val(name), path("${name}.svim.filtered.vcf"), emit: svim_vcf

    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    svim alignment ./ ${bam} ${params.ref} \
        --sample ${name} \
        --min_sv_size 50 \
        --max_sv_size 100000 \
        --types DEL,INS,INV,DUP:TANDEM,DUP:INT \
        --minimum_depth 4  &&

    bcftools sort variants.vcf --temp-dir ./ -o sorted.variants.vcf
    bcftools view -i 'QUAL >= 5' sorted.variants.vcf -o sorted.QUAL.variants.vcf
    bcftools view -i 'FILTER="PASS"' sorted.QUAL.variants.vcf -o ${name}.svim.filtered.vcf
    """
}

/* ------------------------- SvAnna ------------------------- */

//SvAnna - Efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. For installation see: {https://github.com/TheJacksonLaboratory/SvAnna}
//SvAnna annotation doesn't work well with the combiSV output. Seems to ignore all deletions. Better to run SvAnna on each sv caller output independently

process annotate_svim_vcf {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/SvAnna/svim/${name}", mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    tuple val(name), path(vcf_file) 
    path(HPO_terms)

    output:
    path("${name}.svim.SvAnna.*")

    script:
    """
    cp ${params.HPO_terms} phenotype_terms.txt

    source /opt/conda/bin/activate lrs_analysis
    java -jar /opt/svanna-cli-1.0.4/svanna-cli-1.0.4.jar prioritize \
        --data-directory=/opt/svanna-data/ \
        --output-format html,tsv,vcf \
        --out-dir ./ \
        --prefix ${name}.svim.SvAnna \
        --n-threads ${task.cpus} \
        --uncompressed-output \
        --report-top-variants 50 \
        --vcf ${vcf_file} `cat phenotype_terms.txt`
    """
}

/* ------------------------- Sniffles2 ------------------------- */

// To install sniffles2 either do: {conda install sniffles=2.0}

process sniffles2 {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/sniffles/${name}", mode: 'copy', pattern: "${name}.sniffles.filtered.vcf" 
    cpus 8
    memory '16 GB'

    input:
    tuple val(name), file(bam), file(bai)

    output:
    tuple val(name), path("${name}.sniffles.filtered.vcf"), emit: sniffles_vcf

    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    sniffles \
        --threads ${task.cpus} \
        --sample-id ${name} \
        --reference ${params.ref} \
        --minsvlen 50 \
        --symbolic \
        --cluster-merge-pos 150 \
        --tandem-repeats ${tr_bed} \
        --input ${bam} \
        --vcf ${name}.sniffles.vcf &&

    sed '/.:0:0:0:NULL/d' ${name}.sniffles.vcf > tmp.vcf
    mv tmp.vcf ${name}.sniffles.vcf
    bcftools view -i 'FILTER="PASS"' ${name}.sniffles.vcf -o ${name}.sniffles.filtered.vcf
    """
}

/* ------------------------- SvAnna ------------------------- */

//SvAnna - Efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. For installation see: {https://github.com/TheJacksonLaboratory/SvAnna}
//SvAnna annotation doesn't work well with the combiSV output. Seems to ignore all deletions. Better to run SvAnna on each SV caller output independently

process annotate_sniffles2_vcf {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/SvAnna/sniffles/${name}", mode: 'copy' 
    cpus 4
    memory '8 GB'

    input:
    tuple val(name), path(vcf_file)
    path(HPO_terms)

    output:
    path("${name}.sniffles.SvAnna.*")

    script:
    """
    cp ${params.HPO_terms} phenotype_terms.txt

    source /opt/conda/bin/activate lrs_analysis
    java -jar /opt/svanna-cli-1.0.4/svanna-cli-1.0.4.jar prioritize \
        --data-directory=/opt/svanna-data/ \
        --output-format html,tsv,vcf \
        --out-dir ./ \
        --prefix ${name}.sniffles.SvAnna \
        --n-threads ${task.cpus} \
        --uncompressed-output \
        --report-top-variants 50 \
        --vcf ${vcf_file} `cat phenotype_terms.txt`
    """
}

/* ------------------------- combiSV ------------------------- */

//combiSV combines sample SV calls from sniffles, svim, and cuteSV to generate a high-confident call SV VCF per sample. {see: https://github.com/ndierckx/combiSV.git}

process mergeCalls {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/CombiSV/${name}", mode: 'copy' 
    cpus 2
    memory '4 GB'

    input:
    tuple val(name), file(cutesv_vcf), file(svim_vcf), file(sniffles_vcf)

    output:
    tuple val(name), path("${name}.combiSV.merged.consensus.vcf"), emit: consensus_vcf

    script:
    // Use 'name' variable for output file naming
    """
    perl /opt/combiSV/combiSV2.3.pl \
        -sniffles ${sniffles_vcf} \
        -cutesv ${cutesv_vcf} \
        -svim ${svim_vcf} \
        -o ${name}.combiSV.merged.consensus.vcf

    # modify combisv vcf
    awk 'BEGIN{OFS="\t"} /^##/ || /^#/ {print; next} {split(\$8,a,";"); for(i in a) {if(a[i] ~ /^SVTYPE=/) {svtype=substr(a[i],8)}; if(a[i] ~ /^SVLEN=/) {svlen=substr(a[i],7); if(svtype=="DEL") {gsub(svlen,"-"svlen,\$8);}}} \$1="chr"\$1; \$4="N"; \$5="<"svtype">"; print}' ${name}.combiSV.merged.consensus.vcf > output.vcf
    mv output.vcf ${name}.combiSV.merged.consensus.vcf

    """
}

/* ------------------------- AnnotSV ------------------------- */

//AnnotSV - An integrated tool for Structural Variations annotation and ranking. For installation see: {https://github.com/lgmgeo/AnnotSV}
//Annotates the output of consensus calls from combiSV and emits tsv,logs from Annotsv

process annotate_mergedCalls {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/SV_Calls/Annotation/AnnotSV/${name}", mode: 'copy'
    cpus 12
    memory '16 GB'

    input:
    tuple val(name), path(merged_vcf)

    output:
    path("${name}.combiSV.merged.consensus.Annotsv.tsv")
    path("${name}.combiSV.merged.consensus.Annotsv.log") // Also capture the log file

    script:
    """
    source /opt/conda/bin/activate lrs_analysis
    AnnotSV \
        -annotationsDir ${params.annotationsDir} \
        -SVinputFile ${merged_vcf} \
        -outputFile ./${name}.combiSV.merged.consensus.Annotsv.tsv \
        -snvIndelPASS 1 \
        -SVminSize 50 \
        -genomeBuild GRCh38 \
        > ${name}.combiSV.merged.consensus.Annotsv.log
    """
}

/* ------------------------- NanoRepeat ------------------------- */

//NanoRepeat: quantification of Short Tandem Repeats (STRs) from long-read sequencing data
//To install and set up NanoRepeat, see https://github.com/WGLab/NanoRepeat#installation

process nanoRepeat {
    label 'lrs_analysis'
    publishDir "${params.output_dir}/STRs/${name}", mode: 'copy'
    cpus 4
    memory '12 GB'

    input:
    tuple val(name), file(bam), file(bai)

    output:
    path("${name}.details")
    path("${name}.NanoRepeat_output.tsv") // capture summary file

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
        -o ${name}
    """
}

workflow {
    // --- MODIFIED WORKFLOW LOGIC ---
    // Scenario 1: --str only (no --consensus)
    if (params.str && !params.consensus) {
        println("Running STR analysis only (--str specified without --consensus).")
        if (params.STR_regions != "" && params.STR_regions != null) {
            STR_analysis_results = nanoRepeat(input_ch)
        } else {
            println("ERROR: --str flag used but --STR_regions not provided. Skipping STR quantification.")
            // exit 1, "ERROR: --str flag requires --STR_regions parameter."
        }
    }
    // Scenarios 2, 3, and 4 involve SV calling, so group them
    else {
        // --- Always run SV calling pipelines if not in 'STR only' mode ---
        cutesv_calling_results = cutesv(input_ch)
        svim_calling_results = svim(input_ch)
        sniffles_calling_results = sniffles2(input_ch)

        // --- Always run SvAnna annotation if HPO terms are provided (if not 'STR only') ---
        if (params.HPO_terms != "" && params.HPO_terms != null) {
            HPO_terms_ch = Channel.value(file(params.HPO_terms)) // Create channel for HPO terms file
            cutesv_vcf_annotation_results = annotate_cutesv_vcf(cutesv_calling_results, HPO_terms_ch)
            svim_vcf_annotation_results = annotate_svim_vcf(svim_calling_results, HPO_terms_ch)
            sniffles_vcf_annotation_results = annotate_sniffles2_vcf(sniffles_calling_results, HPO_terms_ch)
        } else {
            println("--HPO_terms not provided: Skipping SvAnna annotation of VCFs from individual SV callers [CuteSV, Sniffles2, SVIM]")
        }

        // --- Now differentiate based on --consensus ---
        if (params.consensus) {
            // Scenario 2: --consensus only (no --str)
            // Scenario 3: --consensus AND --str
            println("Running consensus SV analysis (--consensus specified).")

            // --- Merge VCFs ---
            merged_vcfs = cutesv_calling_results.join(svim_calling_results, by: 0).join(sniffles_calling_results, by: 0)
            combisv_results = mergeCalls(merged_vcfs)

            // --- Annotate Merged VCFs if annotationsDir is provided ---
            if (params.annotationsDir != "" && params.annotationsDir != null) {
                 Annotsv_annotation_results = annotate_mergedCalls(combisv_results)
            } else {
                 println("--annotationsDir not provided: Skipping AnnotSV annotation of consensus SVs")
            }

            // --- If --str is ALSO set (Scenario 3), run nanoRepeat ---
            if (params.str) {
                println("Running STR analysis (--str specified).")
                if (params.STR_regions != "" && params.STR_regions != null) {
                    STR_analysis_results = nanoRepeat(input_ch)
                } else {
                     println("ERROR: --str flag used but --STR_regions not provided. Skipping STR quantification.")
                     // exit 1, "ERROR: --str flag requires --STR_regions parameter."
                }
            } else {
                 println("Skipping STR analysis (only --consensus specified, not --str).")
            }

        } else {
            // Scenario 4: Neither --consensus nor --str
            println("Running individual SV callers and optional SvAnna annotation only (neither --consensus nor --str specified).")
            println("Skipping consensus merging (combiSV) and merged annotation (AnnotSV).")
            println("Skipping STR analysis (nanoRepeat).")
        }
    }
}