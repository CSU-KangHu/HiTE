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
    nextflow run main.nf --genome ../demo/genome.fa --thread 48 --out_dir ../demo/test --plant 0

    Mandatory arguments:
      --genome      Genome assembly path (format: fasta, fa, and fna)
      --out_dir      Output directory; It is recommended to use a new directory to avoid automatic deletion of important files.
    General options:
      --work_dir                        Temporary work dir. default = [ /tmp ]
      --threads                         Input thread num. default = [ 10 ]
      --chunk_size                      The chunk size of large genome, default = [ 400 MB ]
      --plant                           Is it a plant genome, 1: true, 0: false. default = [ 1 ]
      --curated_lib                     Provide a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome. We recommend using TE libraries from Repbase. default = [ None ]
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
      --is_denovo_nonltr                Whether to detect non-ltr de novo, 1: true, 0: false. default = [ 0 ]
      --debug                           Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ 0 ]
      --use_HybridLTR                   Whether to use HybridLTR to identify LTRs, 1: true, 0: false. default = [1]
      --use_NeuralTE                    Whether to use NeuralTE to classify TEs, 1: true, 0: false. default = [1]
      --is_wicker                       Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [0]
      --is_output_LTR_lib               Whether to output LTR library. default = [1]

      --flanking_len                    The flanking length of candidates to find the true boundaries, default = [ 50 ]
      --fixed_extend_base_threshold     The length of variation can be tolerated during pairwise alignment, default = [ 1000 ]
      --tandem_region_cutoff            Cutoff of the candidates regarded as tandem region, default = [ 0.5 ]
      --max_repeat_len                  The maximum length of a single repeat, default = [ 30000 ]
      --min_TE_len                      The minimum TE length, default = [ 80 ]
      --chrom_seg_length                The length of genome segments, default = [ 500000 ]
    """.stripIndent()
}


def printSetting() {
    log.info"""
    ====================================Parameter settings========================================
      [Setting] Reference sequences / assemblies path = [ $params.genome ]
      [Setting] work_dir = [ $params.work_dir ]
      [Setting] threads = [ $params.threads ]
      [Setting] Is classified = [ $params.classified ]
      [Setting] Is remove nested TE = [ $params.remove_nested ]
      [Setting] Is getting domain = [ $params.domain ]
      [Setting] The neutral mutation rate (per bp per ya) = [ $params.miu ]
      [Setting] The chunk size of large genome = [ $params.chunk_size ] MB
      [Setting] Is plant genome = [ $params.plant ]
      [Setting] Curated library = [ $params.curated_lib ]
      [Setting] recover = [ $params.recover ]
      [Setting] annotate = [ $params.annotate ]
      [Setting] BM_RM2 = [ $params.BM_RM2 ]
      [Setting] BM_EDTA = [ $params.BM_EDTA ]
      [Setting] BM_HiTE = [ $params.BM_HiTE ]
      [Setting] skip_HiTE = [ $params.skip_HiTE ]
      [Setting] is_denovo_nonltr = [ $params.is_denovo_nonltr ]
      [Setting] debug = [ $params.debug ]
      [Setting] Output Directory = [ $params.out_dir ]
      [Setting] use_HybridLTR = [ $params.use_HybridLTR ]
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

if (!params.out_dir){
    exit 1, "--out_dir option not specified!"
}

printSetting()
projectDir = workflow.projectDir
genome_name = file(params.genome).getName()
out_genome = "${params.out_dir}/${genome_name}"

filePrefix = genome_name.substring(0, genome_name.lastIndexOf('.'))
out_genome_rename = "${params.out_dir}/${filePrefix}.rename.fa"

//parameters of HiTE
tmp_output_dir = file(params.out_dir).toAbsolutePath()
chrom_seg_length = "${params.chrom_seg_length}"
chunk_size = "${params.chunk_size}"
fixed_extend_base_threshold = "${params.fixed_extend_base_threshold}"
max_repeat_len = "${params.max_repeat_len}"
min_TE_len = "${params.min_TE_len}"
flanking_len = "${params.flanking_len}"
tandem_region_cutoff = "${params.tandem_region_cutoff}"
recover = "${params.recover}"
plant = "${params.plant}"
curated_lib = "${params.curated_lib}"
classified = "${params.classified}"
domain = "${params.domain}"
annotate = "${params.annotate}"
debug = "${params.debug}"
threads = "${params.threads}"
is_denovo_nonltr = "${params.is_denovo_nonltr}"
miu = "${params.miu}"
ref = "${params.genome}"
use_HybridLTR = "${params.use_HybridLTR}"
use_NeuralTE = "${params.use_NeuralTE}"
is_output_LTR_lib = "${params.is_output_LTR_lib}"
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
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/coarseBoundary"

    input:
    tuple path(cut_ref), path(prev_TE), path(ref)

    output:
    path "longest_repeats_${ref_index}.flanked.fa"

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
     --recover ${recover} --debug ${debug} -w ${params.work_dir}
    """
}


process TIR {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/TIR"

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
    --prev_TE ${prev_TE} \
    --all_low_copy_tir ${tmp_output_dir}/tir_low_copy.fa \
    -w ${params.work_dir}
    """

}

process Helitron {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/Helitron"

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
    --prev_TE ${prev_TE} --all_low_copy_helitron ${tmp_output_dir}/helitron_low_copy.fa \
    -w ${params.work_dir}
    """
}

process Non_LTR {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/Non_LTR"

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
    --split_ref_dir ${ref_chr} \
    --is_denovo_nonltr ${is_denovo_nonltr} \
    -r ${ref} --prev_TE ${prev_TE} \
    --all_low_copy_non_ltr ${tmp_output_dir}/non_ltr_low_copy.fa \
    -w ${params.work_dir}
    """
}

process GetPrevTEs {
    label 'process_low'

    output:
    path "prev_TE.fa",    emit: ch_pre_tes

    script:
    cores = task.cpus
    """
     touch prev_TE.fa

     if [ -f "${curated_lib}" ]; then
        cat ${curated_lib} >> prev_TE.fa
     fi
    """
}

process GetPrevTEsForLTR {
    label 'process_low'

    input:
    path tir
    path helitron
    path non_ltr
    path other

    output:
    path "tmp_lib.fa",    emit: ch_tmp_lib

    script:
    cores = task.cpus
    """
     if [ -f "${curated_lib}" ]; then
        cat ${curated_lib} > tmp_lib.fa
     fi

     cat ${tir} ${helitron} ${non_ltr} ${other} >> tmp_lib.fa
    """
}

process GenomeClean {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/GenomeClean"

    input:
    path ref

    output:
    path "genome.fa.clean",    emit: ch_genome_clean

    script:
    cores = task.cpus
    """
    genome_clean.py \
     -i ${ref} \
     -o genome.fa.clean \
     -t ${cores} \
     -w ${params.work_dir}
    """
}

process OtherTE {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/OtherTE"

    input:
    path ref

    output:
    path "confident_other.fa",    emit: ch_others

    script:
    cores = task.cpus
    """
    judge_Other_transposons.py \
     -t ${cores} \
     --recover ${recover} -r ${ref} \
     --min_TE_len ${min_TE_len} \
     -w ${params.work_dir}
    """
}

process LTR {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/LTR"

    input:
    path ref
    path tmp_lib

    output:
    path "confident_ltr_cut.fa",    emit: ch_LTRs
    path "intact_LTR.list",    emit: ch_LTR_pass_list
    path "intact_LTR.fa",    emit: ch_intact_LTR
    path "intact_LTR.fa.classified",    emit: ch_intact_classified_LTR
    path "chr_name.map",    emit: chr_name_map

    script:
    cores = task.cpus
    """
    judge_LTR_transposons.py \
     -g ${ref} -t ${cores} --recover ${recover} \
     --use_HybridLTR ${use_HybridLTR} \
     --use_NeuralTE ${use_NeuralTE} --miu ${miu} \
     --is_wicker ${is_wicker} --is_output_lib ${is_output_LTR_lib} \
     -w ${params.work_dir} --prev_TE ${tmp_lib}
    """
}

process UnwrapNested {
    cpus { params.threads ?: 1 }

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
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/BuildLib"

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
     --is_wicker ${is_wicker} -w ${params.work_dir}
    """
}

process annotate_chunk {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/annotate_chunk/${genome_name}"

    input:
    val genome_name
    tuple path(chunk_fasta), path(panTE_lib)

    output:
    tuple val(genome_name), path("${genome_name}_${ref_index}.gff"), path("${genome_name}_${ref_index}.full_length.gff"), emit: ch_annotated_chunks

    script:
    (full, ref_index) = (chunk_fasta =~ /genome.cut(\d+)\.fa/)[0]
    """
    pan_annotate_genome.py --threads ${task.cpus} --panTE_lib ${panTE_lib} \
    --reference ${chunk_fasta} --genome_name ${genome_name}_${ref_index} -w ${params.work_dir}
    """
}

process merge_annotations {
    storeDir "${params.out_dir}/merge_annotations/${genome_name}"

    input:
    tuple val(genome_name), path(annotated_chunks), path(full_length_annotated_chunks)

    output:
    tuple val(genome_name), path("${genome_name}_merged.gff"), path("${genome_name}_merged.full_length.gff"), emit: ch_merged_annotations

    script:
    """
    cat ${annotated_chunks} > ${genome_name}_merged.gff
    cat ${full_length_annotated_chunks} > ${genome_name}_merged.full_length.gff
    """
}

process recover_split_annotation {
    storeDir "${params.out_dir}/recover_split_annotation/${genome_name}"

    input:
    path genome_path
    tuple val(genome_name), path(input_gff), path(full_length_input_gff), path(te_lib)

    output:
    path "${genome_name}.sorted.gff", emit: ch_gff
    path "${genome_name}.sorted.gff.tbl", emit: ch_gff_tbl
    path "${genome_name}.full_length.sorted.gff", emit: ch_full_length_gff
    path "${genome_name}.full_length.sorted.gff.tbl", emit: ch_full_length_gff_tbl

    script:
    """
    recover_split_annotation.py ${input_gff} ${genome_name}.gff
    bedtools sort -i ${genome_name}.gff > ${genome_name}.sorted.gff
    get_summary_count.sh ${te_lib} ${genome_name}.sorted.gff ${genome_path}

    recover_split_annotation.py ${full_length_input_gff} ${genome_name}.full_length.gff
    bedtools sort -i ${genome_name}.full_length.gff > ${genome_name}.full_length.sorted.gff
    get_summary_count.sh ${te_lib} ${genome_name}.full_length.sorted.gff ${genome_path} -gff3
    """
}


process Benchmarking {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/Benchmarking"

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
         -r ${ref} --species ${species} --EDTA_home ${EDTA_home} --recover ${recover} -w ${params.work_dir}
        """
    } else {
        """
        benchmarking.py \
         --BM_RM2 ${BM_RM2} \
         --BM_HiTE ${BM_HiTE} \
         -t ${cores} --TE_lib ${TE_lib} \
         -r ${ref} --species ${species} --recover ${recover} -w ${params.work_dir}
        """
    }
}



/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

workflow {
    out_dir = file(params.out_dir).toAbsolutePath().toString()

    // split genome into chunks
    Channel.fromPath(params.genome, type: 'any', checkIfExists: true).set{ ch_genome }

    if (!params.skip_HiTE) {
        if ( !file(params.genome).exists() )
                exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"

        // Step0: dependency check
        dependencies = Channel.from(params.dependencies)
        EnvCheck(dependencies) | view { "$it" }

        GenomeClean(ch_genome)

        ch_genome = GenomeClean.out.ch_genome_clean

        // Step2: homology Non-LTR identification
        OtherTE(ch_genome)

        // get split genome
        SplitGenome(ch_genome)

        //get identified TEs
        GetPrevTEs()

        // After splitting the genome into chunks, we need to associate each chunk with the pre_tes.
        ch_cut_genomes = SplitGenome.out.cut_genomes.flatten()
        ch_cut_genomes_combined = ch_cut_genomes
            .combine(GetPrevTEs.out.ch_pre_tes)
            .combine(ch_genome)

        // Step3: Coarse Boundary Repeat Sequence Identification
        ch_coarse_TEs = coarseBoundary(ch_cut_genomes_combined)

        // Combine the coarse boundary results with other outputs and input them into the subsequent module.
        ch_coarse_TEs_combined = ch_coarse_TEs
            .combine(GetPrevTEs.out.ch_pre_tes)
            .combine(SplitGenome.out.ref_chr)
            .combine(ch_genome)

        // Step4: TIR identification
        TIR(ch_coarse_TEs_combined)

        // Step5: Helitron identification
        Helitron(ch_coarse_TEs_combined)

        // Step6: non-LTR identification
        Non_LTR(ch_coarse_TEs_combined)

        // Merge all chunks and store them in the output directory.
        all_tirs = TIR.out.ch_TIRs.collectFile(name: "${out_dir}/confident_tir.fa")
        all_helitrons = Helitron.out.ch_Helitrons.collectFile(name: "${out_dir}/confident_helitron.fa")
        all_non_ltrs = Non_LTR.out.ch_Non_LTRs.collectFile(name: "${out_dir}/confident_non_ltr.fa")
        all_others = OtherTE.out.ch_others.collectFile(name: "${out_dir}/confident_other.fa")

        GetPrevTEsForLTR(all_tirs, all_helitrons, all_non_ltrs, all_others)
        // Step1: LTR identification
        LTR(ch_genome, GetPrevTEsForLTR.out.ch_tmp_lib)
        all_ltrs = LTR.out.ch_LTRs.collectFile(name: "${out_dir}/confident_ltr_cut.fa")
        all_LTR_pass_list = LTR.out.ch_LTR_pass_list.collectFile(name: "${out_dir}/intact_LTR.list")
        all_chr_name_map = LTR.out.chr_name_map.collectFile(name: "${out_dir}/chr_name.map")

        // Step7: Build TE library
        BuildLib(all_ltrs, all_tirs, all_helitrons, all_non_ltrs, all_others)

        ch_TEs = BuildLib.out.ch_TEs.collectFile(name: "${out_dir}/confident_TE.cons.fa")
        ch_classified_TE = BuildLib.out.ch_classified_TE.collectFile(name: "${out_dir}/TE_merge_tmp.fa.classified")

        if (params.annotate){
            // Step8: Genome annotation
            annotate_input_channel = ch_cut_genomes.combine(ch_TEs)

            annotate_chunk_out = annotate_chunk("HiTE", annotate_input_channel)
            annotate_chunk_out.ch_annotated_chunks
                .groupTuple(by: 0)
                .set { grouped_chunks }

            merge_annotations_out = merge_annotations(grouped_chunks)
            ch_merge_annotations_out = merge_annotations_out.ch_merged_annotations
            ch_merge_annotations_out.combine(ch_TEs).set { recover_split_annotation_input_channel }
            recover_split_annotation_out = recover_split_annotation(ch_genome, recover_split_annotation_input_channel)
            recover_split_annotation_out.ch_gff.collectFile(name: "${out_dir}/HiTE.sorted.gff")
            recover_split_annotation_out.ch_gff_tbl.collectFile(name: "${out_dir}/HiTE.sorted.gff.tbl")
            recover_split_annotation_out.ch_full_length_gff.collectFile(name: "${out_dir}/HiTE.full_length.sorted.gff")
            recover_split_annotation_out.ch_full_length_gff_tbl.collectFile(name: "${out_dir}/HiTE.full_length.sorted.gff.tbl")
        }
    } else {
        Channel.fromPath("${out_dir}/confident_TE.cons.fa", type: 'any', checkIfExists: true).set{ ch_TEs }
    }

    // Step10: conduct benchmarking
    Benchmarking(ch_TEs, ch_genome)
    Benchmarking.out.ch_BM_RM2.collectFile(name: "${out_dir}/BM_RM2.log")
    Benchmarking.out.ch_BM_HiTE.collectFile(name: "${out_dir}/BM_HiTE.log")
    Benchmarking.out.ch_BM_EDTA.collectFile(name: "${out_dir}/BM_EDTA.log")
}


/*
========================================================================================
    THE END
========================================================================================
*/
