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
ch_EAHelitron = "${projectDir}/bin/EAHelitron-master"
ch_ltrfinder = "${projectDir}/bin/LTR_FINDER_parallel-master"
genome_name = file(params.genome).getName()
out_genome = "${params.outdir}/${genome_name}"

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
    tag "${reference}"

    label 'process_low'

    input:
    path reference
    val ch_module
    val tmp_output_dir
    val chrom_seg_length
    val chunk_size

    output:
    path "${reference}.cut*.fa"

    script:
    cores = task.cpus
    """
    cp ${reference} ${tmp_output_dir}

    python3 ${ch_module}/split_genome_chunks.py \
     -g ${reference} --tmp_output_dir ${tmp_output_dir} --chrom_seg_length ${chrom_seg_length} \
     --chunk_size ${chunk_size}
    """
}

process coarseBoundary {
    tag "${cut_reference}"

    label 'process_high'

    input:
    path cut_reference
    val tmp_output_dir
    val fixed_extend_base_threshold
    val max_repeat_len
    val flanking_len
    val tandem_region_cutoff
    val reference


    output:
    path "longest_repeats_*.flanked.fa"

    script:
    cores = task.cpus
    ref_name = file(reference).getName()
    (full, ref_index) = (cut_reference =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/coarse_boundary.py \
     -g ${cut_reference} --tmp_output_dir ${tmp_output_dir} \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} -r ${reference}

    ## Since nextflow will look for output files in the work directory, we need to copy the script output files to the work directory.
    cp ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa ./
    """
}

process TIR {
    tag "${cut_reference}"

    label 'process_high'

    input:
    path cut_reference
    path longest_repeats_flanked
    val tmp_output_dir
    val plant
    val flanking_len
    val tandem_region_cutoff
    val reference


    output:
    path "confident_tir_*.fa"

    script:
    cores = task.cpus
    ref_name = file(reference).getName()
    (full, ref_index) = (cut_reference =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/judge_TIR_transposons.py \
     -g ${cut_reference} --seqs ${longest_repeats_flanked} \
     -t ${cores} --TRsearch_dir ${tools_module}  \
     --tmp_output_dir ${tmp_output_dir} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} --plant ${plant}

    cp ${tmp_output_dir}/confident_tir_${ref_index}.fa ./
    """
}

process Helitron {
    tag "${cut_reference}"

    label 'process_high'

    input:
    path cut_reference
    path longest_repeats_flanked
    val tmp_output_dir
    val flanking_len
    val EAHelitron
    val reference


    output:
    path "confident_helitron_*.fa"

    script:
    cores = task.cpus
    ref_name = file(reference).getName()
    (full, ref_index) = (cut_reference =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/judge_Helitron_transposons.py \
     -g ${cut_reference} --seqs ${longest_repeats_flanked} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --flanking_len ${flanking_len} --EAHelitron ${EAHelitron} \
     --ref_index ${ref_index}

    cp ${tmp_output_dir}/confident_helitron_${ref_index}.fa ./
    """
}

process OtherTE {
    tag "${cut_reference}"

    label 'process_high'

    input:
    path cut_reference
    path longest_repeats_flanked
    val tmp_output_dir
    val reference


    output:
    path "confident_other_*.fa"

    script:
    cores = task.cpus
    ref_name = file(reference).getName()
    (full, ref_index) = (cut_reference =~ /${ref_name}.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/judge_Other_transposons.py \
     --seqs ${longest_repeats_flanked} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --query_coverage 0.8 --subject_coverage 0 \
     --ref_index ${ref_index}

    cp ${tmp_output_dir}/confident_other_${ref_index}.fa ./
    """
}

process LTR {
    tag "${reference}"

    label 'process_high'

    input:
    path reference
    val ltrfinder_home
    val tmp_output_dir


    output:
    path "confident_ltr_cut.fa"

    script:
    cores = task.cpus
    """
    python3 ${ch_module}/judge_LTR_transposons.py \
     -g ${reference} --ltrfinder_home ${ltrfinder_home} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --recover 0

    cp ${tmp_output_dir}/confident_ltr_cut.fa ./
    """
}

process UnwrapNested {
    tag "${reference}"

    label 'process_high'

    input:
    path reference
    path confident_ltr_cut
    path confident_tir
    path confident_helitron
    path confident_other
    val tmp_output_dir
    val global_flanking_filter
    val remove_nested
    val test_home


    output:
    path "confident_TE.fa"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/remove_nested.py \
     -g ${reference} --confident_ltr_cut ${confident_ltr_cut} \
     --confident_tir ${confident_tir} \
     --confident_helitron ${confident_helitron} \
     --confident_other ${confident_other} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --global_flanking_filter ${global_flanking_filter} \
     --remove_nested ${remove_nested} --test_home ${test_home}

    cp ${tmp_output_dir}/confident_TE.fa ./
    """
}

process BuildLib {
    tag "${confident_TE}"

    label 'process_high'

    input:
    path confident_TE
    val tmp_output_dir
    val classified
    val TEClass_home
    val debug
    val reference

    output:
    path "confident_TE.cons.fa*"


    script:
    cores = task.cpus
    ref_name = file(reference).getName()
    """
    python3 ${ch_module}/get_nonRedundant_lib.py \
     --confident_TE ${confident_TE} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --classified ${classified} --TEClass_home ${TEClass_home} \
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
    Channel.fromPath(params.genome, type: 'any', checkIfExists: true).set{ ch_genome }
    cut_genomes = splitGenome(ch_genome, ch_module, params.outdir, params.chrom_seg_length, params.chunk_size)
    ch_cut_genomes = cut_genomes.flatten()

    //coarse-grained Boundary identification
    longest_repeats = coarseBoundary(ch_cut_genomes, params.outdir,
        params.fixed_extend_base_threshold, params.max_repeat_len,
        params.flanking_len, params.tandem_region_cutoff, out_genome)

    //TIR identification
    ch_tirs = TIR(ch_cut_genomes, longest_repeats, params.outdir, params.plant, params.flanking_len, params.tandem_region_cutoff, out_genome).collectFile(name: "${params.outdir}/confident_tir.fa")

    //Helitron identification
    ch_helitrons = Helitron(ch_cut_genomes, longest_repeats, params.outdir, params.flanking_len, ch_EAHelitron, out_genome).collectFile(name: "${params.outdir}/confident_helitron.fa")

    //Other identification
    ch_others = OtherTE(ch_cut_genomes, longest_repeats, params.outdir, out_genome).collectFile(name: "${params.outdir}/confident_other.fa")
    //test(ch_others) | view { "$it" }

    //LTR identification
    ch_ltrs = LTR(ch_genome, ch_ltrfinder, params.outdir)
    //test(ch_ltrs) | view { "$it" }

    //Unwrap nested TE
    ch_TE = UnwrapNested(ch_genome, ch_ltrs, ch_tirs, ch_helitrons, ch_others, params.outdir, params.global_flanking_filter, params.remove_nested, ch_module)
    //test(ch_TE) | view { "$it" }

    //Build TE library
    ch_lib = BuildLib(ch_TE, params.outdir, params.classified, ch_classification, params.debug, out_genome)
    //test(ch_lib) | view { "$it" }
}


/*
========================================================================================
    THE END
========================================================================================
*/
