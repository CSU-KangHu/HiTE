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
      --recover                         Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ 0 ]
      --miu                             The neutral mutation rate (per bp per ya). default = [ 1.3e-8 ]
      --classified                      Whether to classify TE models, HiTE uses RepeatClassifier from RepeatModeler to classify TEs, 1: true, 0: false. default = [ 1 ]
      --remove_nested                   Whether to remove nested TE, 1: true, 0: false. default = [ 1 ]
      --domain                          Whether to obtain TE domains, HiTE uses RepeatPeps.lib from RepeatMasker to obtain TE domains, 1: true, 0: false. default = [ 0 ]
      --annotate                        Whether to annotate the genome using the TE library generated, 1: true, 0: false. default = [ 0 ]
      --BM_RM2                          Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false. default = [ 0 ]
      --BM_EDTA                         Whether to conduct benchmarking of EDTA, 1: true, 0: false. default = [ 0 ]
      --BM_HiTE                         Whether to conduct benchmarking of HiTE, 1: true, 0: false. default = [ 0 ]
      --EDTA_home                       When conducting benchmarking of EDTA, you will be asked to input EDTA home path.
      --species                         Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).
      --skip_HiTE                       Whether to skip_HiTE, 1: true, 0: false. default = [ 0 ]
      --is_prev_mask                    Whether to mask current genome used the TEs detected in previous iteration, 1: true, 0: false. default = [ 1 ]
      --is_denovo_nonltr                Whether to detect non-ltr de novo, 1: true, 0: false. default = [ 0 ]
      --debug                           Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ 0 ]
      --use_NeuralTE                    Whether to use NeuralTE to classify TEs, 1: true, 0: false. default = [1]
      --is_wicker                       Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [0]

      --flanking_len                    The flanking length of candidates to find the true boundaries, default = [ 50 ]
      --fixed_extend_base_threshold     The length of variation can be tolerated during pairwise alignment, default = [ 1000 ]
      --tandem_region_cutoff            Cutoff of the candidates regarded as tandem region, default = [ 0.5 ]
      --max_repeat_len                  The maximum length of a single repeat, default = [ 30000 ]
      --chrom_seg_length                The length of genome segments, default = [ 500000 ]
    """.stripIndent()
}


def printSetting() {
    log.info"""
    ====================================Parameter settings========================================
      [Setting] Reference sequences / assemblies path = [ $params.genome ]
      [Setting] Is classified = [ $params.classified ]
      [Setting] Is remove nested TE = [ $params.remove_nested ]
      [Setting] Is getting domain = [ $params.domain ]
      [Setting] The neutral mutation rate (per bp per ya) = [ $params.miu ]
      [Setting] The chunk size of large genome = [ $params.chunk_size ] MB
      [Setting] Is plant genome = [ $params.plant ]
      [Setting] recover = [ $params.recover ]
      [Setting] annotate = [ $params.annotate ]
      [Setting] BM_RM2 = [ $params.BM_RM2 ]
      [Setting] BM_EDTA = [ $params.BM_EDTA ]
      [Setting] BM_HiTE = [ $params.BM_HiTE ]
      [Setting] skip_HiTE = [ $params.skip_HiTE ]
      [Setting] is_prev_mask = [ $params.is_prev_mask ]
      [Setting] is_denovo_nonltr = [ $params.is_denovo_nonltr ]
      [Setting] debug = [ $params.debug ]
      [Setting] Output Directory = [ $params.outdir ]
      [Setting] use_NeuralTE = [ $params.use_NeuralTE ]
      [Setting] is_wicker = [ $params.is_wicker ]

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
HSDIR = "${projectDir}/bin/HelitronScanner/TrainingSet"
HSJAR = "${projectDir}/bin/HelitronScanner/HelitronScanner.jar"
sh_dir = "${projectDir}/bin"
member_script_path = "${projectDir}/tools/make_fasta_from_blast.sh"
subset_script_path = "${projectDir}/tools/ready_for_MSA.sh"
lib_module = "${projectDir}/library"
ch_EAHelitron = "${projectDir}/bin/EAHelitron-master"
ch_ltrfinder = "${projectDir}/bin/LTR_FINDER_parallel-master"
ch_ltrharvest = "${projectDir}/bin/LTR_HARVEST_parallel"
ch_NeuralTE = "${projectDir}/bin/NeuralTE"
ch_protein = "${projectDir}/library/RepeatPeps.lib"
genome_name = file(params.genome).getName()
out_genome = "${params.outdir}/${genome_name}"

filePrefix = genome_name.substring(0, genome_name.lastIndexOf('.'))
out_genome_rename = "${params.outdir}/${filePrefix}.rename.fa"

//parameters of HiTE
tmp_output_dir = "${params.outdir}"
chrom_seg_length = "${params.chrom_seg_length}"
chunk_size = "${params.chunk_size}"
fixed_extend_base_threshold = "${params.fixed_extend_base_threshold}"
max_repeat_len = "${params.max_repeat_len}"
flanking_len = "${params.flanking_len}"
tandem_region_cutoff = "${params.tandem_region_cutoff}"
recover = "${params.recover}"
plant = "${params.plant}"
classified = "${params.classified}"
domain = "${params.domain}"
annotate = "${params.annotate}"
is_prev_mask = "${params.is_prev_mask}"
debug = "${params.debug}"
is_denovo_nonltr = "${params.is_denovo_nonltr}"
miu = "${params.miu}"
ref = "${params.genome}"
use_NeuralTE = "${params.use_NeuralTE}"
is_wicker = "${params.is_wicker}"
//parameters of Evaluation
BM_RM2 = "${params.BM_RM2}"
BM_EDTA = "${params.BM_EDTA}"
BM_HiTE = "${params.BM_HiTE}"
rm2_script = "${projectDir}/bin/get_family_summary_paper.sh"
rm2_strict_script = "${projectDir}/bin/get_family_summary_paper_0.99.sh"
EDTA_home = "${params.EDTA_home}"
species = "${params.species}"

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
    path "genome.cut*.fa"
    

    script:
    cores = task.cpus
    ref_name = ref.getName()
    """
    cp ${ref} ${tmp_output_dir}
    python3 ${ch_module}/split_genome_chunks.py \
     -g ${tmp_output_dir}/${ref_name} --tmp_output_dir ${tmp_output_dir} --chrom_seg_length ${chrom_seg_length} \
     --chunk_size ${chunk_size}
    cp ${tmp_output_dir}/genome.cut*.fa ./
    """
}

process coarseBoundary {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(prev_TE)

    output:
    path "longest_repeats_*.flanked.fa"

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/coarse_boundary.py \
     -g ${cut_ref} --tmp_output_dir ${tmp_output_dir} \
     --prev_TE ${tmp_output_dir}/prev_TE.fa \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} -r ${out_genome} --recover ${recover} --is_prev_mask ${is_prev_mask} --debug ${debug}

    ## Since nextflow will look for output files in the work directory, we need to copy the script output files to the work directory.
    cp ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa ./
    """
}

process TE_identification {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(prev_TE)

    output:
    path "confident_tir_*.fa", emit: ch_tirs
    path "confident_helitron_*.fa", emit: ch_helitrons
    path "confident_non_ltr_*.fa", emit: ch_non_ltrs

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]
    """
    # Step1: De novo TE searching
    python3 ${ch_module}/coarse_boundary.py \
     -g ${cut_ref} --tmp_output_dir ${tmp_output_dir} \
     --prev_TE ${tmp_output_dir}/prev_TE.fa \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} -r ${out_genome} --recover ${recover} --is_prev_mask ${is_prev_mask} --debug ${debug}

    ## Since nextflow will look for output files in the work directory, we need to copy the script output files to the work directory.
    cp ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa ./

    # Step2: TIR identification
    python3 ${ch_module}/judge_TIR_transposons.py \
    -g ${cut_ref} --seqs ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa \
    -t ${cores} --TRsearch_dir ${tools_module}  \
    --tmp_output_dir ${tmp_output_dir} \
    --tandem_region_cutoff ${tandem_region_cutoff} \
    --ref_index ${ref_index} \
    --subset_script_path ${subset_script_path} \
    --plant ${plant} \
    --flanking_len ${flanking_len} \
    --recover ${recover} \
    --debug ${debug} \
    -r ${ref} \
    --split_ref_dir ${tmp_output_dir}/ref_chr

    cp ${tmp_output_dir}/confident_tir_${ref_index}.fa ./

    # Step3: Helitron identification
    python3 ${ch_module}/judge_Helitron_transposons.py \
    --seqs ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa -r ${ref} -t ${cores} \
    --tmp_output_dir ${tmp_output_dir} --HSDIR ${HSDIR} --HSJAR ${HSJAR} \
    --sh_dir ${sh_dir} --EAHelitron ${ch_EAHelitron} \
    --subset_script_path ${subset_script_path} \
    --ref_index ${ref_index} --flanking_len ${flanking_len} \
    --recover ${recover} --debug ${debug} --split_ref_dir ${tmp_output_dir}/ref_chr

    cp ${tmp_output_dir}/confident_helitron_${ref_index}.fa ./

    # Step4: non-LTR identification
    python3 ${ch_module}/judge_Non_LTR_transposons.py \
    --seqs ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa -t ${cores} \
    --subset_script_path ${subset_script_path} \
    --tmp_output_dir ${tmp_output_dir} \
    --library_dir ${lib_module} \
    --recover ${recover} \
    --plant ${plant} \
    --debug ${debug} \
    --flanking_len ${flanking_len} \
    --ref_index ${ref_index} \
    --is_denovo_nonltr ${is_denovo_nonltr} \
    -r ${ref}

    cp ${tmp_output_dir}/confident_non_ltr_${ref_index}.fa ./
    """
}

def groupTuple = {
    input_ch.collectFileGroups { file -> file.baseName.replaceAll("[^\\d]", "") }
}

process TIR {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(prev_TE)
    path(lrf)

    output:
    path "confident_tir_*.fa"


    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]

    script:
    """
    python3 ${ch_module}/judge_TIR_transposons.py \
    -g ${cut_ref} --seqs ${lrf} \
    -t ${cores} --TRsearch_dir ${tools_module}  \
    --tmp_output_dir ${tmp_output_dir} \
    --tandem_region_cutoff ${tandem_region_cutoff} \
    --ref_index ${ref_index} \
    --subset_script_path ${subset_script_path} \
    --plant ${plant} \
    --flanking_len ${flanking_len} \
    --recover ${recover} \
    --debug ${debug} \
    -r ${ref} \
    --split_ref_dir ${tmp_output_dir}/ref_chr

    cp ${tmp_output_dir}/confident_tir_${ref_index}.fa ./
    """

}

process Helitron {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(lrf)

    output:
    path "confident_helitron_*.fa"

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]

    script:
    """
    python3 ${ch_module}/judge_Helitron_transposons.py \
    --seqs ${lrf} -r ${ref} -t ${cores} \
    --tmp_output_dir ${tmp_output_dir} --HSDIR ${HSDIR} --HSJAR ${HSJAR} \
    --sh_dir ${sh_dir} --EAHelitron ${ch_EAHelitron} \
    --subset_script_path ${subset_script_path} \
    --ref_index ${ref_index} --flanking_len ${flanking_len} \
    --recover ${recover} --debug ${debug} --split_ref_dir ${tmp_output_dir}/ref_chr

    cp ${tmp_output_dir}/confident_helitron_${ref_index}.fa ./
    """
}

process Non_LTR {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(lrf)

    output:
    path "confident_non_ltr_*.fa"

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]

    script:
    """
    python3 ${ch_module}/judge_Non_LTR_transposons.py \
    --seqs ${lrf} -t ${cores} \
    --subset_script_path ${subset_script_path} \
    --tmp_output_dir ${tmp_output_dir} \
    --library_dir ${lib_module} \
    --recover ${recover} \
    --plant ${plant} \
    --debug ${debug} \
    --flanking_len ${flanking_len} \
    --ref_index ${ref_index} \
    -r ${ref}

    cp ${tmp_output_dir}/confident_non_ltr_${ref_index}.fa ./
    """
}


process merge_ltr_other {
    label 'process_low'

    input:
    path ltr
    path other

    output:
    path "prev_TE.fa"

    script:
    cores = task.cpus
    """
     cat ${ltr} > prev_TE.fa
     cat ${other} >> prev_TE.fa
     cp prev_TE.fa ${tmp_output_dir}/prev_TE.fa
    """
}

process OtherTE {
    tag "${ref}"

    label 'process_high'

    input:
    path ref

    output:
    path "confident_other.fa"

    script:
    cores = task.cpus
    """
    python3 ${ch_module}/judge_Other_transposons.py \
     -r ${ref} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --library_dir ${lib_module} --recover ${recover}

    cp ${tmp_output_dir}/confident_other.fa ./
    """
}

process LTR {
    tag "${ref}"

    label 'process_high'

    input:
    path ref

    output:
    path "confident_ltr_cut.fa.cons"

    script:
    cores = task.cpus
    """
    python3 ${ch_module}/judge_LTR_transposons.py \
     -g ${ref} --ltrharvest_home ${ch_ltrharvest} --ltrfinder_home ${ch_ltrfinder} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --recover ${recover} --miu ${miu} --use_NeuralTE ${use_NeuralTE} --is_wicker ${is_wicker}\
     --NeuralTE_home ${ch_NeuralTE} --TEClass_home ${ch_classification}

    cp ${tmp_output_dir}/confident_ltr_cut.fa.cons ./
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
    label 'process_single'

    input:
    path ltr
    path tir
    path helitron
    path non_ltr
    path other

    output:
    path "confident_TE.cons.fa"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/get_nonRedundant_lib.py \
     --confident_ltr_cut ${ltr} \
     --confident_tir ${tir} \
     --confident_helitron ${helitron} \
     --confident_non_ltr ${non_ltr} \
     --confident_other ${other} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --test_home ${ch_module} --use_NeuralTE ${use_NeuralTE} --is_wicker ${is_wicker} \
     --NeuralTE_home ${ch_NeuralTE} --TEClass_home ${ch_classification} \
     --domain ${domain} --protein_path ${ch_protein}

    cp ${tmp_output_dir}/confident_TE.cons.fa ./
    """
}

process ClassifyLib {
    label 'process_high'

    input:
    path lib

    output:
    path "${lib}.classified"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/get_classified_lib.py \
     --confident_TE_consensus ${lib} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --classified ${classified} --domain ${domain} --TEClass_home ${ch_classification} \
     --protein_path ${ch_protein} \
     --debug ${debug}
    """
}

process annotate_genome {
    label 'process_high'

    input:
    path lib
    path ref

    output:
    file "output.txt"

    script:
    cores = task.cpus
    """
    python3 ${ch_module}/annotate_genome.py \
     -t ${cores} --classified_TE_consensus ${lib} \
     --tmp_output_dir ${tmp_output_dir} \
     --annotate ${annotate} \
      -r ${ref}

    echo "annotate_genome_finished" > output.txt
    """
}

process benchmarking {
    label 'process_high'

    input:
    path TE_lib
    path ref
    path annotate_dout

    output:
    file "output.txt"

    script:
    cores = task.cpus
    if (file("${EDTA_home}/lib-test.pl").exists()) {
        """
        python3 ${ch_module}/benchmarking.py \
         --tmp_output_dir ${tmp_output_dir} \
         --BM_RM2 ${BM_RM2} --BM_EDTA ${BM_EDTA} --BM_HiTE ${BM_HiTE} \
         -t ${cores} --lib_module ${lib_module} --TE_lib ${TE_lib} \
         --rm2_script ${rm2_script} --rm2_strict_script ${rm2_strict_script} \
         -r ${ref} --species ${species} --EDTA_home ${EDTA_home}

        echo "benchmarking_finished" > output.txt
        """
    } else {
        """
        python3 ${ch_module}/benchmarking.py \
         --tmp_output_dir ${tmp_output_dir} \
         --BM_RM2 ${BM_RM2} \
         --BM_HiTE ${BM_HiTE} \
         -t ${cores} --lib_module ${lib_module} --TE_lib ${TE_lib} \
         --rm2_script ${rm2_script} --rm2_strict_script ${rm2_strict_script} \
         -r ${ref} --species ${species}

        echo "benchmarking_finished" > output.txt
        """
    }
}

process CleanLib {
    label 'process_low'

    input:
    path lib

    output:
    stdout


    script:
    """
    python3 ${ch_module}/clean_lib.py \
     --tmp_output_dir ${tmp_output_dir} \
     --debug ${debug}
    """
}

process BM_RM2 {
    tag "${TE}"

    label 'process_high'

    input:
    path TE
    path curatedLib
    path rm2_script

    output:
    path "${TE}.out"
    path "res.log"


    script:
    cores = task.cpus
    """
    RepeatMasker -lib ${curatedLib} -nolow -pa ${cores} ${TE}
    mkdir rm2_test
    cd rm2_test && rm -rf * && sh ../${rm2_script} ../${TE}.out > ../res.log
    """
}

process BM_EDTA {
    tag "${TE}"

    label 'process_high'

    input:
    path TE
    path curatedLib
    val reference
    val EDTA_home

    output:
    path "repbase.out"
    path "HiTE.out"
    path "HiTE.out.*"


    script:
    cores = task.cpus
    """
    RepeatMasker -e ncbi -pa ${cores} -q -no_is -norna -nolow -div 40 -lib ${curatedLib} -cutoff 225 ${tmp_output_dir}/genome.rename.fa
    mv ${tmp_output_dir}/genome.rename.fa.out repbase.out
    cp repbase.out ${tmp_output_dir}/

    RepeatMasker -e ncbi -pa ${cores} -q -no_is -norna -nolow -div 40 -lib ${TE} -cutoff 225 ${tmp_output_dir}/genome.rename.fa
    mv ${tmp_output_dir}/genome.rename.fa.out HiTE.out
    cp HiTE.out ${tmp_output_dir}/

    perl ${EDTA_home}/lib-test.pl -genome ${tmp_output_dir}/genome.rename.fa -std repbase.out -tst HiTE.out -cat Total
    """
}

// allow multi-thread
process test {
    tag "${cut_reference}"

    label 'process_low'

    input:
    val input1

    output:
    stdout

    script:
    cores = task.cpus


    """
    # genome.fa.cut0.fa
    echo "${input1}"
    """

}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

workflow {
    if (!params.skip_HiTE) {
        if ( !file(params.genome).exists() )
                exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"

        // dependency check
        dependencies = Channel.from(params.dependencies)
        EnvCheck(dependencies) | view { "$it" }

        // split genome into chunks
        Channel.fromPath(params.genome, type: 'any', checkIfExists: true).set{ ch_genome }

        //LTR identification
        ch_ltrs = LTR(ch_genome)

        //homology Non-LTR identification
        ch_others = OtherTE(ch_genome)

        //get identified TEs
        ch_pre_tes = merge_ltr_other(ch_ltrs, ch_others)

        // get split genome
        cut_genomes = splitGenome(ch_genome)
        ch_cut_genomes = cut_genomes.flatten()
        ch_cut_genomes_combined = ch_cut_genomes.combine(ch_pre_tes)

        // TIR, Helitron, non-LTR identification
        TE_identification(ch_cut_genomes_combined)
        //test(ch_ouput) | view { "$it" }
        ch_tirs = TE_identification.out.ch_tirs
        ch_helitrons = TE_identification.out.ch_helitrons
        ch_non_ltrs = TE_identification.out.ch_non_ltrs

        all_tirs = ch_tirs.collectFile(name: "${params.outdir}/confident_tir.fa")
        all_helitrons = ch_helitrons.collectFile(name: "${params.outdir}/confident_helitron.fa")
        all_non_ltrs = ch_non_ltrs.collectFile(name: "${params.outdir}/confident_non_ltr.fa")

        //Build TE library
        ch_lib = BuildLib(ch_ltrs, all_tirs, all_helitrons, all_non_ltrs, ch_others)

//         //Classify TE library
//         ch_lib.splitFasta(by: params.classify_chunk_size, file:true).set { ch_fasta }
//         ch_classified_lib = ClassifyLib(ch_fasta)
//         ch_final = ch_classified_lib.collectFile(name: "${params.outdir}/confident_TE.cons.fa.classified")
        //test(ch_lib) | view { "$it" }
        annotate_out = annotate_genome(ch_lib, ch_genome)

    }

    if (params.skip_HiTE)
        Channel.fromPath("${params.outdir}/confident_TE.cons.fa", type: 'any', checkIfExists: true).set{ ch_lib }
    bm_out = benchmarking(ch_lib, ch_genome, annotate_out)
    CleanLib(bm_out)

}


/*
========================================================================================
    THE END
========================================================================================
*/
