#!/usr/bin/env nextflow
/*
========================================================================================
    HiTE
========================================================================================
    Github : https://github.com/CSU-KangHu/HiTE
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

def echo(text) {
  println text
}

def helpMessage() {
    log.info"""
    HiTE - Nextflow PIPELINE (v$workflow.manifest.version)
    =================================
    Usage:
    The typical command is as follows:
    nextflow run main.nf --genome ../demo/genome.fa --thread 48 --outdir ../demo/test --plant 0

    Mandatory arguments:
      --genome      Genome assembly path (format: fasta, fa, and fna)
      --outdir      Output directory; It is recommended to use a new directory to avoid automatic deletion of important files.
    General options:
      --chunk_size                      The chunk size of large genome, default = [ 400 MB ]
      --plant                           Is it a plant genome, 1: true, 0: false. default = [ 1 ]
      --remove_nested                   Whether to unwrap the nested TE, 1: true, 0: false. default = [ 1 ]
      --classified                      Whether to classify TE models, HiTE uses RepeatClassifier from RepeatModeler to classify TEs, 1: true, 0: false. default = [ 1 ]
      --recover                         Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ 0 ]
      --debug                           Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ 0 ]
      --flanking_len                    The flanking length of candidates to find the true boundaries, default = [ 50 ]
      --fixed_extend_base_threshold     The length of variation can be tolerated during pairwise alignment, default = [ 1000 ]
      --tandem_region_cutoff            Cutoff of the candidates regarded as tandem region, default = [ 0.5 ]
      --max_repeat_len                  The maximum length of a single repeat, default = [ 30000 ]
      --chrom_seg_length                The length of genome segments, default = [ 500000 ]
      --global_flanking_filter          Whether to filter false positives by global flanking alignment, significantly reduce false positives but require more memory, especially when inputting a large genome. 1: true (require more memory), 0: false. default = [ 1 ]
    """.stripIndent()
}

def printSetting() {
    log.info"""
    ====================================Parameter settings========================================
      [Setting] Reference sequences / assemblies path = [ $params.genome ]
      [Setting] The chunk size of large genome = [ $params.chunk_size ] MB
      [Setting] Is plant genome = [ $params.plant ]
      [Setting] Remove nested = [ $params.remove_nested ]
      [Setting] Global flanking filter = [ $params.global_flanking_filter ]
      [Setting] recover = [ $params.recover ]
      [Setting] debug = [ $params.debug ]
      [Setting] Output Directory = [ $params.outdir ]

      [Setting] Fixed extend bases threshold = [ $params.fixed_extend_base_threshold ]
      [Setting] Flanking length of TE = [ $params.flanking_len ]
      [Setting] Cutoff of the repeat regarded as tandem sequence = [ $params.tandem_region_cutoff ]
      [Setting] Maximum length of TE = [ $params.max_repeat_len ]
      [Setting] The length of genome segments = [ $params.chrom_seg_length ]
    """.stripIndent()
}

echo "Welcome to HiTE nextflow!"

// check input =====================================================================
// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (!params.genome){
    exit 1, "--genome option not specified!"
}

if (!params.outdir){
    exit 1, "--outdir option not specified!"
}

printSetting()
projectDir = workflow.projectDir
ch_module = "${projectDir}/module"
ch_classification = "${projectDir}/classification"
tools_module = "${projectDir}/tools"
lib_module = "${projectDir}/library"
ch_EAHelitron = "${projectDir}/bin/EAHelitron-master"
ch_ltrfinder = "${projectDir}/bin/LTR_FINDER_parallel-master"
genome_name = file(params.genome).getName()
out_genome = "${params.outdir}/${genome_name}"

//parameters of HiTE
tmp_output_dir = "${params.outdir}"
chrom_seg_length = "${params.chrom_seg_length}"
chunk_size = "${params.chunk_size}"
fixed_extend_base_threshold = "${params.fixed_extend_base_threshold}"
max_repeat_len = "${params.max_repeat_len}"
flanking_len = "${params.flanking_len}"
tandem_region_cutoff = "${params.tandem_region_cutoff}"
plant = "${params.plant}"
global_flanking_filter = "${params.global_flanking_filter}"
remove_nested = "${params.remove_nested}"
classified = "${params.classified}"
debug = "${params.debug}"

// Check all tools work well
process EnvCheck {
    tag "envcheck"
    errorStrategy 'terminate'

    label 'process_low'

    input:
    val dependency

    output:
    stdout

    script:
    """
    echo "---Check if ${dependency} is installed."
    command -v ${dependency} >/dev/null 2>&1 || { echo >&2 "${dependency} is required but not installed.  Aborting."; exit 1; }
    """
}

// Check all tools work well
process splitGenome {
    tag "${ref}"

    label 'process_low'

    input:
    path ref

    output:
    path "${ref}.cut*.fa"

    script:
    cores = task.cpus
    """
    cp ${ref} ${tmp_output_dir}

    python3 ${ch_module}/split_genome_chunks.py \
     -g ${ref} --tmp_output_dir ${tmp_output_dir} --chrom_seg_length ${chrom_seg_length} \
     --chunk_size ${chunk_size}
    """
}

process coarseBoundary {
    tag "${cut_ref}"

    label 'process_high'

    input:
    path cut_ref

    output:
    path "longest_repeats_*.flanked.fa"

    script:
    cores = task.cpus
    ref_name = file(out_genome).getName()
    (full, ref_index) = (cut_ref =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/coarse_boundary.py \
     -g ${cut_ref} --tmp_output_dir ${tmp_output_dir} \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} -r ${out_genome}

    ## Since nextflow will look for output files in the work directory, we need to copy the script output files to the work directory.
    cp ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa ./
    """
}

process TIR {
    tag "${cut_ref}"

    label 'process_high'

    input:
    path cut_ref
    path lrf


    output:
    path "confident_tir_*.fa"

    script:
    cores = task.cpus
    ref_name = file(out_genome).getName()
    (full, ref_index) = (cut_ref =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/judge_TIR_transposons.py \
     -g ${cut_ref} --seqs ${lrf} \
     -t ${cores} --TRsearch_dir ${tools_module}  \
     --tmp_output_dir ${tmp_output_dir} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} --plant ${plant}

    cp ${tmp_output_dir}/confident_tir_${ref_index}.fa ./
    """
}

process Helitron {
    tag "${cut_ref}"

    label 'process_high'

    input:
    path cut_ref
    path lrf


    output:
    path "confident_helitron_*.fa"

    script:
    cores = task.cpus
    ref_name = file(out_genome).getName()
    (full, ref_index) = (cut_ref =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/judge_Helitron_transposons.py \
     -g ${cut_ref} --seqs ${lrf} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --flanking_len ${flanking_len} --EAHelitron ${ch_EAHelitron} \
     --ref_index ${ref_index}

    cp ${tmp_output_dir}/confident_helitron_${ref_index}.fa ./
    """
}

process OtherTE {
    tag "${cut_ref}"

    label 'process_high'

    input:
    path cut_ref
    path lrf


    output:
    path "confident_other_*.fa"

    script:
    cores = task.cpus
    ref_name = file(out_genome).getName()
    (full, ref_index) = (cut_ref =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/judge_Other_transposons.py \
     --seqs ${lrf} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --query_coverage 0.8 --subject_coverage 0 \
     --ref_index ${ref_index} --library_dir ${lib_module}

    cp ${tmp_output_dir}/confident_other_${ref_index}.fa ./
    """
}

process LTR {
    tag "${ref}"

    label 'process_high'

    input:
    path ref

    output:
    path "confident_ltr_cut.fa"

    script:
    cores = task.cpus
    """
    python3 ${ch_module}/judge_LTR_transposons.py \
     -g ${ref} --ltrfinder_home ${ch_ltrfinder} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --recover 0

    cp ${tmp_output_dir}/confident_ltr_cut.fa ./
    """
}

process UnwrapNested {
    tag "${ref}"

    label 'process_high'

    input:
    path ref
    path ltr
    path tir
    path helitron
    path other


    output:
    path "confident_TE.fa"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/remove_nested.py \
     -g ${ref} --confident_ltr_cut ${ltr} \
     --confident_tir ${tir} \
     --confident_helitron ${helitron} \
     --confident_other ${other} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --global_flanking_filter ${global_flanking_filter} \
     --remove_nested ${remove_nested} --test_home ${ch_module}

    cp ${tmp_output_dir}/confident_TE.fa ./
    """
}

process BuildLib {
    tag "${TE}"

    label 'process_high'

    input:
    path TE

    output:
    path "confident_TE.cons.fa*"


    script:
    cores = task.cpus
    ref_name = file(out_genome).getName()
    """
    python3 ${ch_module}/get_nonRedundant_lib.py \
     --confident_TE ${TE} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --classified ${classified} --TEClass_home ${ch_classification} \
     --debug ${debug} --ref_name {ref_name}

    cp ${tmp_output_dir}/confident_TE.cons.fa* ./
    """
}


// allow multi-thread
process test {
    tag "${cut_reference}"

    label 'process_low'

    input:
    val cut_reference

    output:
    stdout

    script:
    cores = task.cpus


    """
    # genome.fa.cut0.fa
    echo "${cut_reference}"
    """

}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

workflow {

    if ( !file(params.genome).exists() )
        exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"


    // dependency check
    dependencies = Channel.from(params.dependencies)
    EnvCheck(dependencies) | view { "$it" }

    // split genome into chunks
    Channel.fromPath(params.genome, type: 'any', checkIfExists: true).set{ ch_g }
    cut_genomes = splitGenome(ch_g)
    ch_cut_g = cut_genomes.flatten()

    //coarse-grained Boundary identification
    longest_repeats = coarseBoundary(ch_cut_g)

    //TIR identification
    ch_tirs = TIR(ch_cut_g, longest_repeats).collectFile(name: "${params.outdir}/confident_tir.fa")

    //Helitron identification
    ch_h = Helitron(ch_cut_g, longest_repeats).collectFile(name: "${params.outdir}/confident_helitron.fa")

    //Other identification
    ch_o = OtherTE(ch_cut_g, longest_repeats).collectFile(name: "${params.outdir}/confident_other.fa")
    //test(ch_o) | view { "$it" }

    //LTR identification
    ch_ltrs = LTR(ch_g)
    //test(ch_ltrs) | view { "$it" }

    //Unwrap nested TE
    ch_TE = UnwrapNested(ch_g, ch_ltrs, ch_tirs, ch_h, ch_o)
    //test(ch_TE) | view { "$it" }

    //Build TE library
    ch_lib = BuildLib(ch_TE)
    //test(ch_lib) | view { "$it" }
}


/*
========================================================================================
    THE END
========================================================================================
*/
