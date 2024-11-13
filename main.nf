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
      --curated_lib                     Provide a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome. We recommend using TE libraries from Repbase. default = [ None ]
      --recover                         Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ 0 ]
      --intact_anno                     Whether to generate annotation of full-length TEs, 1: true, 0: false. default = [ 0 ]
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
      [Setting] Curated library = [ $params.curated_lib ]
      [Setting] recover = [ $params.recover ]
      [Setting] annotate = [ $params.annotate ]
      [Setting] intact_anno = [ $params.intact_anno ]
      [Setting] BM_RM2 = [ $params.BM_RM2 ]
      [Setting] BM_EDTA = [ $params.BM_EDTA ]
      [Setting] BM_HiTE = [ $params.BM_HiTE ]
      [Setting] skip_HiTE = [ $params.skip_HiTE ]
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
curated_lib = "${params.curated_lib}"
classified = "${params.classified}"
domain = "${params.domain}"
annotate = "${params.annotate}"
debug = "${params.debug}"
is_denovo_nonltr = "${params.is_denovo_nonltr}"
miu = "${params.miu}"
ref = "${params.genome}"
use_NeuralTE = "${params.use_NeuralTE}"
is_wicker = "${params.is_wicker}"
search_struct = "${params.search_struct}"
//parameters of Evaluation
BM_RM2 = "${params.BM_RM2}"
BM_EDTA = "${params.BM_EDTA}"
BM_HiTE = "${params.BM_HiTE}"
EDTA_home = "${params.EDTA_home}"
species = "${params.species}"

// 检查依赖工具的Process
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
process SplitGenome {
    tag "${ref}"

    label 'process_low'

    input:
    path ref

    output:
    path "genome.cut*.fa",    emit: cut_genomes
    path "ref_chr",    emit: ref_chr

    script:
    cores = task.cpus
    """
    split_genome_chunks.py \
     -g ${ref} --chrom_seg_length ${chrom_seg_length} --chunk_size ${chunk_size}
    """
}

process coarseBoundary {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(prev_TE), path(ref)

    output:
    path "longest_repeats_*.flanked.fa"

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]
    """
    coarse_boundary.py \
     -g ${cut_ref} --prev_TE ${prev_TE} \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --ref_index ${ref_index} --flanking_len ${flanking_len} \
     --tandem_region_cutoff ${tandem_region_cutoff} -r ${ref} \
     --recover ${recover} --debug ${debug}
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
    coarse_boundary.py \
     -g ${cut_ref} --tmp_output_dir ${tmp_output_dir} \
     --prev_TE ${tmp_output_dir}/prev_TE.fa \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} -r ${out_genome} --recover ${recover} --debug ${debug}

    ## Since nextflow will look for output files in the work directory, we need to copy the script output files to the work directory.
    cp ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa ./

    # Step2: TIR identification
    judge_TIR_transposons.py \
    -g ${cut_ref} --seqs ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa \
    -t ${cores} \
    --tmp_output_dir ${tmp_output_dir} \
    --tandem_region_cutoff ${tandem_region_cutoff} \
    --ref_index ${ref_index} \
    --plant ${plant} \
    --flanking_len ${flanking_len} \
    --recover ${recover} \
    --debug ${debug} \
    -r ${ref} \
    --split_ref_dir ${tmp_output_dir}/ref_chr

    cp ${tmp_output_dir}/confident_tir_${ref_index}.fa ./

    # Step3: Helitron identification
    judge_Helitron_transposons.py \
    --seqs ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa -r ${ref} -t ${cores} \
    --tmp_output_dir ${tmp_output_dir} \
    --ref_index ${ref_index} --flanking_len ${flanking_len} \
    --recover ${recover} --debug ${debug} --split_ref_dir ${tmp_output_dir}/ref_chr

    cp ${tmp_output_dir}/confident_helitron_${ref_index}.fa ./

    # Step4: non-LTR identification
    judge_Non_LTR_transposons.py \
    --seqs ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa -t ${cores} \
    --tmp_output_dir ${tmp_output_dir} \
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
    tag "${lrf}"

    label 'process_high'

    input:
    tuple path(lrf), path(prev_TE), path(ref_chr), path(ref)


    output:
    path "confident_tir_*.fa",    emit: ch_TIRs


    script:
    cores = task.cpus
    (full, ref_index) = (lrf =~ /longest_repeats_(\d+)\.flanked\.fa/)[0]
    script:
    """
    judge_TIR_transposons.py \
    --seqs ${lrf} -t ${cores} \
    --flanking_len ${flanking_len} --plant ${plant} \
    --tandem_region_cutoff ${tandem_region_cutoff} \
    --ref_index ${ref_index} \
    --recover ${recover} \
    --debug ${debug} \
    -r ${ref} \
    --split_ref_dir ${ref_chr} \
    --prev_TE ${prev_TE}
    """

}

process Helitron {
    tag "${lrf}"

    label 'process_high'

    input:
    tuple path(lrf), path(prev_TE), path(ref_chr), path(ref)

    output:
    path "confident_helitron_*.fa",    emit: ch_Helitrons

    script:
    cores = task.cpus
    (full, ref_index) = (lrf =~ /longest_repeats_(\d+)\.flanked\.fa/)[0]

    script:
    """
    judge_Helitron_transposons.py \
    --seqs ${lrf} -r ${ref} -t ${cores} \
    --ref_index ${ref_index} --flanking_len ${flanking_len} \
    --recover ${recover} --debug ${debug} --split_ref_dir ${ref_chr} \
    --prev_TE ${prev_TE}
    """
}

process Non_LTR {
    tag "${lrf}"

    label 'process_high'

    input:
    tuple path(lrf), path(prev_TE), path(ref_chr), path(ref)

    output:
    path "confident_non_ltr_*.fa",    emit: ch_Non_LTRs


    script:
    cores = task.cpus
    (full, ref_index) = (lrf =~ /longest_repeats_(\d+)\.flanked\.fa/)[0]

    script:
    """
    judge_Non_LTR_transposons.py \
    --seqs ${lrf} -t ${cores} \
    --recover ${recover} \
    --plant ${plant} \
    --debug ${debug} \
    --flanking_len ${flanking_len} \
    --ref_index ${ref_index} \
    --is_denovo_nonltr ${is_denovo_nonltr} \
    -r ${ref} --prev_TE ${prev_TE}
    """
}


process MergeLTROther {
    label 'process_low'

    input:
    path ltr
    path other

    output:
    path "prev_TE.fa",    emit: ch_pre_tes

    script:
    cores = task.cpus
    """
     cat ${ltr} > prev_TE.fa

     if [ -f "${curated_lib}" ]; then
        cat ${curated_lib} >> prev_TE.fa
     fi
    """
}

process OtherTE {
    tag "${ref}"

    label 'process_high'

    input:
    path ref

    output:
    path "confident_other.fa",    emit: ch_others

    script:
    cores = task.cpus
    """
    judge_Other_transposons.py \
     -t ${cores} \
     --recover ${recover} -r ${ref}
    """
}

process LTR {
    tag "${ref}"

    label 'process_high'

    input:
    path ref

    output:
    path "confident_ltr_cut.fa",    emit: ch_LTRs, optional: true
    path "genome.rename.fa.pass.list",    emit: ch_LTR_pass_list, optional: true
    path "chr_name.map",    emit: chr_name_map, optional: true

    script:
    cores = task.cpus
    """
    judge_LTR_transposons.py \
     -g ${ref} \
     -t ${cores} --recover ${recover} --use_NeuralTE ${use_NeuralTE} --miu ${miu} \
     --is_wicker ${is_wicker}
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
    remove_nested.py \
     -g ${ref} --confident_ltr_cut ${ltr} \
     --confident_tir ${tir} \
     --confident_helitron ${helitron} \
     --confident_other ${other} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --global_flanking_filter ${global_flanking_filter} \
     --remove_nested ${remove_nested}

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
    path "confident_TE.cons.fa",    emit: ch_TEs
    path "TE_merge_tmp.fa.classified",    emit: ch_classified_TE


    script:
    cores = task.cpus
    """
    get_nonRedundant_lib.py \
     -t ${cores} \
     --confident_ltr_cut ${ltr} \
     --confident_tir ${tir} \
     --confident_helitron ${helitron} \
     --confident_non_ltr ${non_ltr} \
     --confident_other ${other} \
     --use_NeuralTE ${use_NeuralTE} \
     --domain ${domain} --curated_lib ${curated_lib} \
     --is_wicker ${is_wicker}
    """
}

process IntactTEAnnotation {
    label 'process_high'

    input:
    path ch_TEs
    path all_LTR_pass_list
    path all_chr_name_map
    path all_tirs
    path all_helitrons
    path all_non_ltrs
    path all_others
    path ch_genome
    path ch_classified_TE

    output:
    path "HiTE_intact.sorted.gff3", emit:ch_intact_gff


    script:
    cores = task.cpus
    """
    get_full_length_annotation.py \
     -t ${cores} --ltr_list ${all_LTR_pass_list} \
     --tir_lib ${all_tirs}  \
     --helitron_lib ${all_helitrons} \
     --nonltr_lib ${all_non_ltrs} \
     --other_lib ${all_others} \
     --chr_name_map ${all_chr_name_map} \
     -r ${ch_genome} \
     --search_struct ${search_struct} \
     --classified_TE_path ${ch_classified_TE}
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
    get_classified_lib.py \
     --confident_TE_consensus ${lib} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --classified ${classified} --domain ${domain} \
     --debug ${debug}
    """
}

process AnnotateGenome {
    label 'process_high'

    input:
    path lib
    path ref

    output:
    path "HiTE.out",    emit: ch_HiTE_out, optional: true
    path "HiTE.tbl",    emit: ch_HiTE_tbl, optional: true
    path "HiTE.gff",    emit: ch_HiTE_gff, optional: true

    script:
    cores = task.cpus
    """
    annotate_genome.py \
     -t ${cores} --classified_TE_consensus ${lib} \
     --annotate ${annotate} \
      -r ${ref}
    """
}

process Benchmarking {
    label 'process_high'

    input:
    path TE_lib
    path ref

    output:
    path "BM_RM2.log",    emit: ch_BM_RM2, optional: true
    path "BM_HiTE.log",    emit: ch_BM_HiTE, optional: true
    path "BM_EDTA.log",    emit: ch_BM_EDTA, optional: true

    script:
    cores = task.cpus
    if (file("${EDTA_home}/lib-test.pl").exists()) {
        """
        benchmarking.py \
         --BM_RM2 ${BM_RM2} --BM_EDTA ${BM_EDTA} --BM_HiTE ${BM_HiTE} \
         -t ${cores} --TE_lib ${TE_lib} \
         -r ${ref} --species ${species} --EDTA_home ${EDTA_home}
        """
    } else {
        """
        benchmarking.py \
         --BM_RM2 ${BM_RM2} \
         --BM_HiTE ${BM_HiTE} \
         -t ${cores} --TE_lib ${TE_lib} \
         -r ${ref} --species ${species}
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
    clean_lib.py \
     --tmp_output_dir ${tmp_output_dir} \
     --debug ${debug}
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
    outdir = file(params.outdir).toAbsolutePath().toString()

    // split genome into chunks
    Channel.fromPath(params.genome, type: 'any', checkIfExists: true).set{ ch_genome }

    if (!params.skip_HiTE) {
        if ( !file(params.genome).exists() )
                exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"

        // Step0: dependency check
        dependencies = Channel.from(params.dependencies)
        EnvCheck(dependencies) | view { "$it" }

        // Step1: LTR identification
        LTR(ch_genome)

        // Step2: homology Non-LTR identification
        OtherTE(ch_genome)

        //get identified TEs
        MergeLTROther(LTR.out.ch_LTRs, OtherTE.out.ch_others)

        // get split genome
        SplitGenome(ch_genome)

        // After splitting the genome into chunks, we need to associate each chunk with the pre_tes.
        ch_cut_genomes = SplitGenome.out.cut_genomes.flatten()
        ch_cut_genomes_combined = ch_cut_genomes.combine(MergeLTROther.out.ch_pre_tes)
        ch_cut_genomes_combined = ch_cut_genomes_combined.combine(ch_genome)

        // Step3: Coarse Boundary Repeat Sequence Identification
        ch_coarse_TEs = coarseBoundary(ch_cut_genomes_combined)

        // Combine the coarse boundary results with other outputs and input them into the subsequent module.
        ch_coarse_TEs_combined = ch_coarse_TEs.combine(MergeLTROther.out.ch_pre_tes)
        ch_coarse_TEs_combined = ch_coarse_TEs_combined.combine(SplitGenome.out.ref_chr)
        ch_coarse_TEs_combined = ch_coarse_TEs_combined.combine(ch_genome)

        // Step4: TIR identification
        TIR(ch_coarse_TEs_combined)

        // Step5: Helitron identification
        Helitron(ch_coarse_TEs_combined)

        // Step6: non-LTR identification
        Non_LTR(ch_coarse_TEs_combined)
        // test(Non_LTR.out.ch_Non_LTRs) | view { "$it" }

        // Merge all chunks and store them in the output directory.
        all_ltrs = LTR.out.ch_LTRs.collectFile(name: "${outdir}/confident_ltr_cut.fa")
        all_tirs = TIR.out.ch_TIRs.collectFile(name: "${outdir}/confident_tir.fa")
        all_helitrons = Helitron.out.ch_Helitrons.collectFile(name: "${outdir}/confident_helitron.fa")
        all_non_ltrs = Non_LTR.out.ch_Non_LTRs.collectFile(name: "${outdir}/confident_non_ltr.fa")
        all_others = OtherTE.out.ch_others.collectFile(name: "${outdir}/confident_other.fa")
        all_LTR_pass_list = LTR.out.ch_LTR_pass_list.collectFile(name: "${outdir}/genome.rename.fa.pass.list")
        all_chr_name_map = LTR.out.chr_name_map.collectFile(name: "${outdir}/chr_name.map")

        // Step7: Build TE library
        BuildLib(all_ltrs, all_tirs, all_helitrons, all_non_ltrs, all_others)

        ch_TEs = BuildLib.out.ch_TEs.collectFile(name: "${outdir}/confident_TE.cons.fa")
        BuildLib.out.ch_classified_TE.collectFile(name: "${outdir}/TE_merge_tmp.fa.classified")

        // Step8: Genome annotation
        AnnotateGenome(ch_TEs, ch_genome)
        AnnotateGenome.out.ch_HiTE_out.collectFile(name: "${outdir}/HiTE.out")
        AnnotateGenome.out.ch_HiTE_tbl.collectFile(name: "${outdir}/HiTE.tbl")
        AnnotateGenome.out.ch_HiTE_gff.collectFile(name: "${outdir}/HiTE.gff")
    } else {
        Channel.fromPath("${outdir}/confident_TE.cons.fa", type: 'any', checkIfExists: true).set{ ch_TEs }
    }

    if (params.intact_anno){
        // Step9: get full-length TE annotation
        Channel.fromPath("${outdir}/genome.rename.fa.pass.list", type: 'any', checkIfExists: false).set{ all_LTR_pass_list }
        Channel.fromPath("${outdir}/chr_name.map", type: 'any', checkIfExists: false).set{ all_chr_name_map }
        Channel.fromPath("${outdir}/confident_tir.fa", type: 'any', checkIfExists: false).set{ all_tirs }
        Channel.fromPath("${outdir}/confident_helitron.fa", type: 'any', checkIfExists: false).set{ all_helitrons }
        Channel.fromPath("${outdir}/confident_non_ltr.fa", type: 'any', checkIfExists: false).set{ all_non_ltrs }
        Channel.fromPath("${outdir}/confident_other.fa", type: 'any', checkIfExists: false).set{ all_others }
        Channel.fromPath("${outdir}/TE_merge_tmp.fa.classified", type: 'any', checkIfExists: false).set{ ch_classified_TE }
        IntactTEAnnotation(ch_TEs, all_LTR_pass_list, all_chr_name_map, all_tirs, all_helitrons, all_non_ltrs, all_others, ch_genome, ch_classified_TE)
        IntactTEAnnotation.out.ch_intact_gff.collectFile(name: "${outdir}/HiTE_intact.sorted.gff3")
    }

    // Step10: conduct benchmarking
    Benchmarking(ch_TEs, ch_genome)
    Benchmarking.out.ch_BM_RM2.collectFile(name: "${outdir}/BM_RM2.log")
    Benchmarking.out.ch_BM_HiTE.collectFile(name: "${outdir}/BM_HiTE.log")
    Benchmarking.out.ch_BM_EDTA.collectFile(name: "${outdir}/BM_EDTA.log")
}


/*
========================================================================================
    THE END
========================================================================================
*/
