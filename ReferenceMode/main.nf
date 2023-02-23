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
      --thread                          Thread number, default = [the maximum number of cores supported on your machine]
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
      [Setting] Threads = [ $params.thread ]
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

// Check all tools work well
process EnvCheck {
    tag "envcheck"
    errorStrategy 'terminate'

    output:
    stdout

    script:
    """
    echo "### 检查依赖是否安装"
    command -v python >/dev/null 2>&1 || { echo >&2 "python is required but not installed.  Aborting."; exit 1; }
    command -v gt >/dev/null 2>&1 || { echo >&2 "gt is required but not installed.  Aborting."; exit 1; }
    command -v LTR_retriever >/dev/null 2>&1 || { echo >&2 "LTR_retriever is required but not installed.  Aborting."; exit 1; }
    command -v RepeatMasker >/dev/null 2>&1 || { echo >&2 "RepeatMasker is required but not installed.  Aborting."; exit 1; }
    command -v rmblastn >/dev/null 2>&1 || { echo >&2 "rmblastn is required but not installed.  Aborting."; exit 1; }
    command -v RepeatModeler >/dev/null 2>&1 || { echo >&2 "RepeatModeler is required but not installed.  Aborting."; exit 1; }
    command -v bedtools >/dev/null 2>&1 || { echo >&2 "bedtools is required but not installed.  Aborting."; exit 1; }
    echo "### 所有的依赖已经满足"
    """
}

// allow multi-thread
process multi_process {
    label 'process_high'

    output:
    stdout

    script:
    cores = task.cpus
    """
    date; hostname; pwd
    echo "### Check env"
    echo "cpus=${cores}"
    echo "### Check env DONE"
    """
}


workflow {
    //EnvCheck | view { "$it" }
    multi_process() | view { "$it" }
}