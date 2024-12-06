#!/usr/bin/env nextflow
/*
========================================================================================
    panHiTE
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
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

// 定义用户输入参数
params.pan_genomes_dir = ''
params.genome_list = ''
params.genes_dir = '/dev/gene'
params.RNA_dir = '/dev/RNA'
params.out_dir = './output'
params.te_type = 'all'
params.skip_analyze = 0
params.softcore_threshold = 0.8
params.debug = 0
params.threads = 10
params.miu = 1.3e-8
params.all_te_types = ['ltr', 'tir', 'helitron', 'non-ltr', 'all']

// 验证 TE 类型是否合法
if (!params.all_te_types.contains(params.te_type)) {
    error "Invalid TE type: ${params.te_type}. Please choose from ${params.all_te_types}"
}


def helpMessage() {
    log.info"""
    panHiTE - Nextflow PIPELINE (v$workflow.manifest.version)
    =================================
    Usage:
    The typical command is as follows:
    nextflow run panHiTE.nf --pan_genomes_dir xxx --genome_list xxx --genes_dir xxx --RNA_dir xxx --out_dir xxx --threads 40 --skip_analyze 0 --miu 7e-9

    Mandatory arguments:
      --pan_genomes_dir      A directory containing the pan-genomes
      --genome_list          A text file with genome and gene names. Each line represents a pair of genome and gene names, separated by a tab (\t). The genome name is mandatory, while the gene name is optional. If a gene name is provided, the genes_dir parameter must also be specified.
      --out_dir              Output directory
    General options:
      --softcore_threshold   occurrence of core_TE = num_of_genomes, softcore_threshold * num_of_genomes <= softcore_TE < num_of_genomes, 2 <= dispensable_TE < softcore_threshold * num_of_genomes, private_TE = 1. default = [ 0.8 ]
      --genes_dir            A directory containing the gene annotation files, gff format.
      --RNA_dir              A directory containing the RNA-seq files.
      --te_type              Retrieve specific type of TE output [ltr|tir|helitron|non-ltr|all]. default = [ all ]
      --threads              Input thread num. default = [ 10 ]
      --skip_analyze         Whether to skip analyze, only generate panTE library. default = [ 0 ]
      --miu                  The neutral mutation rate (per bp per ya). default = [ 1.3e-8 ]
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

process preprocess_genomes {
    storeDir "${params.out_dir}/preprocess_genomes"

    input:
    path genome_list
    path genes_dir
    path RNA_dir
    path pan_genomes_dir

    output:
    path "genome_metadata.json"

    script:
    """
    echo "Running preprocessing"
    pan_preprocess_genomes.py --genome_list ${genome_list} --genes_dir ${genes_dir} \
    --RNA_dir ${RNA_dir} --pan_genomes_dir ${pan_genomes_dir}
    """
}

// Step 3: HiTE 并行处理每个基因组
process run_hite_single {
    storeDir "${params.out_dir}/run_hite_single/${genome_name}"

    input:
    tuple val(genome_name), val(raw_name), path(reference), val(threads), val(te_type), val(miu), val(debug)

    output:
    tuple val(genome_name), path("intact_LTR.list"), emit: ch_intact_ltr_list, optional: true
    path "confident_ltr.terminal.fa", emit: ch_ltr_terminal, optional: true
    path "confident_ltr.internal.fa", emit: ch_ltr_internal, optional: true
    path "confident_helitron.fa", emit: ch_helitron, optional: true
    path "confident_non_ltr.fa", emit: ch_non_ltr, optional: true
    path "confident_other.fa", emit: ch_other, optional: true
    path "confident_tir.fa", emit: ch_tir, optional: true
    path "confident_TE.cons.fa", emit: ch_te

    script:
    """
    pan_run_hite_single.py --genome_name ${genome_name} --reference ${reference} --threads ${threads} \
    --te_type ${te_type} --miu ${miu} --debug ${debug}
    """
}

// 合并 pan_te.tmp.fa 到 pan_internal.tmp.fa
process merge_terminal_te {
    input:
    path terminal
    path te

    output:
    path "pan_terminal_te.merge.fa"

    script:
    """
    cat ${terminal} ${te} > pan_terminal_te.merge.fa
    """
}

// Step 5: 去冗余
process pan_remove_redundancy {
    input:
    file terminal_tmp
    file internal_tmp
    val threads

    output:
    path "panTE.fa"

    script:
    """
    pan_remove_redundancy.py --pan_terminal_tmp_lib ${terminal_tmp} --pan_internal_tmp_lib ${internal_tmp} --threads ${threads}
    """
}

// Step 3: 注释基因组
process annotate_genomes {
    storeDir "${params.out_dir}/annotate_genomes/${genome_name}"

    input:
    tuple val(genome_name), path(reference), val(threads), path(panTE_lib)

    output:
    tuple val(genome_name), path("${genome_name}.gff"), path("${genome_name}.full_length.gff"), path("${genome_name}.full_length.copies"), emit: annotate_out

    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.gff"

    script:
    """
    pan_annotate_genome.py --threads ${threads} --panTE_lib ${panTE_lib} --reference ${reference} \
    --genome_name ${genome_name}
    """
}

// Step 4: 汇总 TE 数据
process summarize_tes {
    storeDir "${params.out_dir}/summarize_tes"

    input:
    path genome_info_json
    path pan_genomes_dir
    path panTE_lib
    val softcore_threshold

    output:
    path "TE_summary.pdf", emit: ch_te_summary
    path "panHiTE.CorePan_fitmodel.pdf", emit: ch_corepan_model, optional: true
    path "panHiTE.CorePan_fitsmooth.pdf", emit: ch_corepan_smodel, optional: true

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.pdf"

    script:
    """
    pan_summary_TEs.py --genome_info_json ${genome_info_json} --pan_genomes_dir ${pan_genomes_dir} \
    --panTE_lib ${panTE_lib} --softcore_threshold ${softcore_threshold}
    """
}


process pan_gene_te_relation {
    storeDir "${params.out_dir}/pan_gene_te_relation"

    input:
    path genome_info_json

    output:
    path "gene_te_associations.tsv"

    publishDir "${params.out_dir}", mode: 'copy', pattern: "gene_te_associations.tsv"

    script:
    """
    pan_gene_te_relation.py --genome_info_json ${genome_info_json}
    """
}

process pan_detect_de_genes {
    storeDir "${params.out_dir}/pan_detect_de_genes"

    input:
    path genome_info_json
    path gene_te_associations
    val RNA_dir
    val threads

    output:
    path "DE_genes_from_TEs.tsv", emit: ch_de_genes, optional: true
    path "all_gene_TEs_details.tsv", emit: ch_all_genes, optional: true

    publishDir "${params.out_dir}", mode: 'copy', pattern: "DE_genes_from_TEs.tsv, all_gene_TEs_details.tsv"

    script:
    """
    pan_detect_de_genes.py --genome_info_json ${genome_info_json} --gene_te_associations ${gene_te_associations} \
    --RNA_dir ${RNA_dir} --threads ${threads}
    """
}


// 定义工作流
workflow {
    // Step 1: 预处理，生成json格式的输入文件路径
    genome_metadata_out = preprocess_genomes(params.genome_list, params.genes_dir, params.RNA_dir, params.pan_genomes_dir)

    // Step 2: 解析Step1的json文件，准备HiTE输入channel
    genome_info_list = genome_metadata_out
        .flatMap { json_file ->
            def genome_data = file(json_file).text  // 读取文件内容为字符串
            def genome_json = new JsonSlurper().parseText(genome_data) // 使用 JsonSlurper 解析 JSON
            genome_json.genome_info.collect { [it.genome_name, it.raw_name, it.reference, it.gene_gtf, it.RNA_seq] }
        }
    genome_info_list.map { genome_name, raw_name, reference, gene_gtf, RNA_seq ->
            [genome_name, raw_name, reference, params.threads, params.te_type, params.miu, params.debug]
        }.set { hite_input_channel }


    // Step 3: HiTE 并行处理每个基因组
    hite_out = run_hite_single(hite_input_channel)
    // 将每个 Channel 的输出文件收集并合并
    all_terminal = hite_out.ch_ltr_terminal.collectFile(name: "${params.out_dir}/pan_terminal.tmp.fa")
    all_internal = hite_out.ch_ltr_internal.collectFile(name: "${params.out_dir}/pan_internal.tmp.fa")
    all_te = hite_out.ch_te.collectFile(name: "${params.out_dir}/pan_te.tmp.fa")
    intact_ltr_list_channel = hite_out.ch_intact_ltr_list

    // Step 4: 合并 HiTE 的 其他TE和LTR terminal library
    merged_terminal = merge_terminal_te(all_terminal, all_te)
    all_terminal = merged_terminal.collectFile(name: "${params.out_dir}/pan_terminal.tmp.fa")

    // Step 5: 对LTR terminal 和 internal 去冗余，生成panTE library
    panTE_lib = pan_remove_redundancy(all_terminal, all_internal, params.threads)
    panTE_lib = panTE_lib.collectFile(name: "${params.out_dir}/panTE.fa")

    // 准备panTE library和其他参数，作为channel
    annotate_input = genome_info_list.map { genome_name, raw_name, reference, gene_gtf, RNA_seq ->
        [genome_name, reference, params.threads]
    }.combine(panTE_lib).set { annotate_input_channel }

    // Step 6: 并行注释每个基因组
    annotate_out = annotate_genomes(annotate_input_channel)

    // 将注释结果合并到 json 文件中，便于后续统一分析
    genome_info_list.join(annotate_out).map { genome_info ->
        [
            genome_info[0],
            genome_info[1],
            genome_info[2],
            genome_info[3],
            genome_info[4],
            genome_info[5].toString(),
            genome_info[6].toString(),
            genome_info[7].toString()
        ]
    }.join(intact_ltr_list_channel).map { genome_info ->
        [
            genome_name: genome_info[0],
            raw_name   : genome_info[1],
            reference  : genome_info[2],
            gene_gtf   : genome_info[3],
            RNA_seq    : genome_info[4],
            TE_gff        : genome_info[5],
            full_length_TE_gff: genome_info[6],
            full_length_copies: genome_info[7],
            intact_LTR_list: genome_info[8].toString()
        ]
    }.collect().map { data ->
        def jsonContent = "[\n" + data.collect { JsonOutput.toJson(it) }.join(",\n") + "\n]"
        def filePath = "${params.out_dir}/genome_info.json"
        new File(filePath).write(jsonContent)
        return filePath
    }.set { genome_info_json }

    if (!params.skip_analyze) {
         // Step 7: 对检测到的 TE 进行统计分析
        summarize_tes(genome_info_json, params.pan_genomes_dir, panTE_lib, params.softcore_threshold)

        // Step 8: 基因和 TE 关系
        gene_te_associations_out = pan_gene_te_relation(genome_info_json)

        // Step 9: 检测差异表达基因
        detected_de_genes_out = pan_detect_de_genes(genome_info_json, gene_te_associations_out, params.RNA_dir, params.threads)
    }
}
