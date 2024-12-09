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
params.debug = 0
params.threads = 10


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


process pan_generate_bam_for_RNA_seq {
    storeDir "${params.out_dir}/pan_generate_bam_for_RNA_seq"

    input:
    tuple val(genome_name), val(period), path(reference), val(RNA_dir), val(RNA_seq), val(threads)

    output:
    tuple val(genome_name), path("${genome_name}.${period}.output.sorted.bam"), emit: bam_out

    script:
    """
    pan_generate_bam_for_RNA_seq.py --genome_name ${genome_name} --reference ${reference} \
    --RNA_seq '${RNA_seq}' --RNA_dir ${RNA_dir} --threads ${threads}
    """
}



// 定义工作流
workflow {
    // 定义 JSON 文件路径
    genome_metadata_out = Channel.value('/public/home/hpc194701009/ath_pan_genome/pan_genome/watermelon/panLTR_8_watermelon/G42_periods.json')

    // Step 1: 解析输入 JSON
    genome_metadata_out
    .flatMap { json_file ->
        def genome_data = file(json_file).text  // 读取文件内容为字符串
        def genome_json = new JsonSlurper().parseText(genome_data) // 使用 JsonSlurper 解析 JSON
        genome_json.RNA_seq.collect { rna_seq -> // 遍历 RNA_seq 数据
            [
                genome_json.genome_name,
                rna_seq.Status,
                genome_json.reference,
                params.RNA_dir,
                rna_seq,
                params.threads
            ]
        }
    }.set { generate_bam_input_channel }

    //Step 2: 为RNA_seq生成比对bam
    bam_out = pan_generate_bam_for_RNA_seq(generate_bam_input_channel)
//
//
//     // Step 3: 检测差异表达基因
//     // 将 bam 结果合并到 json 文件中, 为pan_detect_de_genes生成输入channel
//     genome_info_list.join(bam_out).map { genome_info ->
//         [
//             genome_name: genome_info[0],
//             reference: genome_info[1],
//             raw_name: genome_info[2],
//             gene_gtf: genome_info[3],
//             RNA_seq    : genome_info[4],
//             bam: genome_info[5].toString()
//         ]
//     }.collect().map { data ->
//         def jsonContent = "[\n" + data.collect { JsonOutput.toJson(it) }.join(",\n") + "\n]"
//         def filePath = "${params.out_dir}/genome_info_for_bam.json"
//         new File(filePath).write(jsonContent)
//         return filePath
//     }.set { genome_info_for_bam_json }

    //pan_detect_de_genes(genome_info_for_bam_json, gene_te_associations_out, params.RNA_dir, params.threads)

}
