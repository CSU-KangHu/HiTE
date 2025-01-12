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
params.min_sample_threshold = 5
params.min_diff_threshold = 10

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
      --min_sample_threshold When identifying differentially expressed genes caused by TE insertions, for the groups with TE insertions and without TE insertions, the number of samples in at least one group must be greater than the min_sample_threshold.
      --min_diff_threshold   The absolute difference in gene expression values between groups with and without TE insertion. diff = min(group1) - max(group2).
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

process pan_preprocess_genomes {
    cpus = 2

    storeDir "${params.out_dir}/pan_preprocess_genomes"

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
    --RNA_dir ${RNA_dir} --pan_genomes_dir ${pan_genomes_dir} > preprocess_genomes.log 2>&1
    """
}

// Step 3: HiTE 并行处理每个基因组
process pan_run_hite_single {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_run_hite_single/${genome_name}"

    input:
    tuple val(genome_name), val(raw_name), path(reference), val(te_type), val(miu), val(debug)

    output:
    tuple val(genome_name), path("intact_LTR.list"), emit: ch_intact_ltr_list
    path "confident_ltr.terminal.fa", emit: ch_ltr_terminal
    path "confident_ltr.internal.fa", emit: ch_ltr_internal
    path "confident_ltr_cut.fa", emit: ch_ltr_cut
    path "intact_LTR.fa", emit: ch_intact_ltr
    path "intact_LTR.fa.classified", emit: ch_intact_classified_ltr
    path "confident_helitron.fa", emit: ch_helitron
    path "confident_non_ltr.fa", emit: ch_non_ltr
    path "confident_other.fa", emit: ch_other
    path "confident_tir.fa", emit: ch_tir
    path "confident_TE.cons.fa", emit: ch_te
    tuple val(genome_name), path("tir_low_copy.fa"), path("helitron_low_copy.fa"), path("non_ltr_low_copy.fa"), emit: ch_low_copy_files

    script:
    cores = task.cpus
    """
    pan_run_hite_single.py --genome_name ${genome_name} --reference ${reference} --threads ${cores} \
    --te_type ${te_type} --miu ${miu} --debug ${debug} > ${genome_name}.run_hite_single.log 2>&1
    """
}

// 合并 pan_te.tmp.fa 到 pan_internal.tmp.fa
process merge_terminal_te {
    cpus = 1

    storeDir "${params.out_dir}/merge_terminal_te"

    input:
    path terminal
    path tir
    path helitron
    path other
    path non_ltr

    output:
    path "pan_terminal_te.merge.fa"

    script:
    """
    cat ${terminal} ${tir} ${helitron} ${other} ${non_ltr} > pan_terminal_te.merge.fa
    """
}

// Step 5: 去冗余
process pan_remove_redundancy {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_remove_redundancy"

    input:
    file merge_te_file

    output:
    path "panTE.fa"

    script:
    cores = task.cpus
    """
    pan_remove_redundancy.py --merge_te_file ${merge_te_file} --threads ${cores} > pan_remove_redundancy.log 2>&1
    """
}

process pan_recover_low_copy_TEs{
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_recover_low_copy_TEs/${genome_name}"

    input:
    tuple val(genome_name), path(tir_low_copy), path(helitron_low_copy), path(non_ltr_low_copy), path(panTE_lib)

    output:
    path "panTE.recover.fa.classified", emit: ch_recover_TEs
    path "get_copies", emit: chr_copies
    path "real_*.fa.cons", emit: chr_raw_real

    script:
    cores = task.cpus
    """
    pan_recover_low_copy_TEs.py --genome_name ${genome_name} --threads ${cores} --tir_low_copy ${tir_low_copy} \
    --helitron_low_copy ${helitron_low_copy} --non_ltr_low_copy ${non_ltr_low_copy} \
    --panTE_lib ${panTE_lib} --genome_list ${params.genome_list} \
    --pan_genomes_dir ${params.pan_genomes_dir} > pan_recover_low_copy_TEs.log 2>&1
    """
}

process pan_merge_TE_recover {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_merge_TE_recover"

    input:
    path panTE_lib
    path panTE_recover_lib

    output:
    path "panTE.merge_recover.fa", emit: ch_panTE_merge
    path "panTE.merge_recover.fa.clstr", emit: ch_clstr

    script:
    cores = task.cpus
    """
    cat ${panTE_lib} ${panTE_recover_lib} > panTE.merge_recover.redundant.fa
    cd-hit-est -aS 0.95 -aL 0.95 -c 0.8 -d 0 -G 0 -g 1 -A 80 -i panTE.merge_recover.redundant.fa \
    -o panTE.merge_recover.fa -T ${cores} -M 0 > pan_merge_TE_recover.log 2>&1
    """
}

// Step 3: 注释基因组
process pan_annotate_genomes {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_annotate_genomes/${genome_name}"

    input:
    tuple val(genome_name), path(reference), path(panTE_lib)

    output:
    tuple val(genome_name), path("${genome_name}.gff"), path("${genome_name}.tbl"), path("${genome_name}.out"), path("${genome_name}.full_length.gff"), path("${genome_name}.full_length.copies"), emit: annotate_out

    script:
    cores = task.cpus
    """
    pan_annotate_genome.py --threads ${cores} --panTE_lib ${panTE_lib} --reference ${reference} \
    --genome_name ${genome_name} > ${genome_name}.annotate_genomes.log 2>&1
    """
}

// Step 4: 汇总 TE 数据
process pan_summarize_tes {
    cpus = 2

    storeDir "${params.out_dir}/pan_summarize_tes"

    input:
    path genome_info_json
    path pan_genomes_dir
    path panTE_lib
    val softcore_threshold

    output:
    path "TE_summary.pdf", emit: ch_te_summary
    path "panHiTE.CorePan_fitmodel.pdf", emit: ch_corepan_model, optional: true
    path "panHiTE.CorePan_fitsmooth.pdf", emit: ch_corepan_smodel, optional: true
    path "panHiTE_fl.CorePan_fitmodel.pdf", emit: ch_corepan_fl_model, optional: true
    path "panHiTE_fl.CorePan_fitsmooth.pdf", emit: ch_corepan_fl_smodel, optional: true
    path "*.txt", emit: ch_corepan_txt, optional: true
    path "TE_PAV.tsv", emit: ch_pav
    path "full_length_TE_PAV.tsv", emit: ch_full_length_pav
    path "Full_length_TEs_Ratio.json", emit: ch_fl_te_ratio
    path "TEs_Ratio.json", emit: ch_te_ratio
    path "Full_length_TE_Coverage.json", emit: ch_fl_te_coverage
    path "TE_Coverage.json", emit: ch_te_coverage
    path "Full_length_TE_Classes_Ratio.json", emit: ch_fl_te_classes_ratio
    path "TE_Classes_Ratio.json", emit: ch_te_classes_ratio
    path "Full_length_TE_Classes_Coverage.json", emit: ch_fl_te_classes_coverage
    path "TE_Classes_Coverage.json", emit: ch_te_classes_coverage
    path "intact_LTR_insert_time.csv", emit: ch_ltr_insert_time

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.pdf"

    script:
    """
    pan_summary_TEs.py --genome_info_json ${genome_info_json} --pan_genomes_dir ${pan_genomes_dir} \
    --panTE_lib ${panTE_lib} --softcore_threshold ${softcore_threshold} > summarize_tes.log 2>&1
    """
}


process pan_gene_te_relation {
    cpus = 2

    storeDir "${params.out_dir}/pan_gene_te_relation"

    input:
    path genome_info_json

    output:
    path "gene_te_associations.tsv"

    publishDir "${params.out_dir}", mode: 'copy', pattern: "gene_te_associations.tsv"

    script:
    """
    pan_gene_te_relation.py --genome_info_json ${genome_info_json} > pan_gene_te_relation.log 2>&1
    """
}


process pan_generate_bam_for_RNA_seq {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_generate_bam_for_RNA_seq"

    input:
    tuple val(genome_name), path(reference), val(RNA_dir), val(RNA_seq)

    output:
    tuple val(genome_name), path("${genome_name}.output.sorted.bam"), emit: bam_out

    script:
    cores = task.cpus
    """
    pan_generate_bam_for_RNA_seq.py --genome_name ${genome_name} --reference ${reference} \
    --RNA_seq '${RNA_seq}' --RNA_dir ${RNA_dir} --threads ${cores}  > ${genome_name}.pan_generate_bam_for_RNA_seq.log 2>&1
    """
}

process pan_detect_de_genes {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/pan_detect_de_genes"

    input:
    path genome_info_for_bam_json
    path gene_te_associations
    val RNA_dir

    output:
    path "DE_genes_from_TEs.tsv", emit: ch_de_genes, optional: true
    path "all_gene_TEs_details.tsv", emit: ch_all_genes, optional: true
    path "DE_genes_from_TEs.pdf", emit: ch_de_genes_pdf, optional: true
    path "TE_express.table", emit: ch_te_express, optional: true
    path "gene_express.table", emit: ch_gene_express, optional: true

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tsv"

    script:
    cores = task.cpus
    """
    pan_detect_de_genes.py --genome_info_for_bam_json ${genome_info_for_bam_json} --gene_te_associations ${gene_te_associations} \
    --min_sample_threshold ${params.min_sample_threshold} --min_diff_threshold ${params.min_diff_threshold} --RNA_dir ${RNA_dir} \
    --threads ${cores} > pan_detect_de_genes.log 2>&1
    """
}


// 定义工作流
workflow {
    // Step 1: 预处理，生成json格式的输入文件路径
    genome_metadata_out = pan_preprocess_genomes(params.genome_list, params.genes_dir, params.RNA_dir, params.pan_genomes_dir)

    // Step 2: 解析Step1的json文件，准备HiTE输入channel
    genome_info_list = genome_metadata_out
        .flatMap { json_file ->
            def genome_data = file(json_file).text  // 读取文件内容为字符串
            def genome_json = new JsonSlurper().parseText(genome_data) // 使用 JsonSlurper 解析 JSON
            genome_json.genome_info.collect { [it.genome_name, it.raw_name, it.reference, it.gene_gtf, it.RNA_seq] }
        }
    genome_info_list.map { genome_name, raw_name, reference, gene_gtf, RNA_seq ->
            [genome_name, raw_name, reference, params.te_type, params.miu, params.debug]
        }.set { hite_input_channel }


    // Step 3: HiTE 并行处理每个基因组
    hite_out = pan_run_hite_single(hite_input_channel)
    // 将每个 Channel 的输出文件收集并合并
    all_te = hite_out.ch_te.collectFile(name: "${params.out_dir}/pan_te.tmp.fa")
    intact_ltr_list_channel = hite_out.ch_intact_ltr_list
    // 将所有的低拷贝TE收集并合并
    low_copy_files_channel = hite_out.ch_low_copy_files

    // Step 4: 对LTR terminal 和 internal 去冗余，生成panTE library
    panTE_lib = pan_remove_redundancy(all_te)
    panTE_lib = panTE_lib.collectFile(name: "${params.out_dir}/panTE.fa")

    // 准备panTE library和其他参数，作为channel
    recover_input = low_copy_files_channel.combine(panTE_lib).set { recover_input_channel }

    // Step 5: 将泛基因组的低拷贝TE 和 panTE 进行聚类，保留和 panTE 序列不一样的低拷贝TE。
    // 从泛基因组中获取拷贝，判断这些低拷贝TE是否是真实TE
    recover_out = pan_recover_low_copy_TEs(recover_input_channel)
    panTE_recover_lib = recover_out.ch_recover_TEs
    panTE_recover_lib = panTE_recover_lib.collectFile(name: "${params.out_dir}/panTE.recover.fa")

    // Step 6: 将恢复的低拷贝 TEs 与 panTE lib 合并
    merge_out = pan_merge_TE_recover(panTE_lib, panTE_recover_lib)
    panTE_merge_lib = merge_out.ch_panTE_merge
    panTE_merge_lib = panTE_merge_lib.collectFile(name: "${params.out_dir}/panTE.merge_recover.fa")

    // 准备panTE library和其他参数，作为channel
    annotate_input = genome_info_list.map { genome_name, raw_name, reference, gene_gtf, RNA_seq ->
        [genome_name, reference]
    }.combine(panTE_merge_lib).set { annotate_input_channel }
 
    // Step 5: 并行注释每个基因组
    annotate_out = pan_annotate_genomes(annotate_input_channel)

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
            genome_info[7].toString(),
            genome_info[8].toString(),
            genome_info[9].toString()
        ]
    }.join(intact_ltr_list_channel).map { genome_info ->
        [
            genome_name: genome_info[0],
            raw_name   : genome_info[1],
            reference  : genome_info[2],
            gene_gtf   : genome_info[3],
            RNA_seq    : genome_info[4],
            TE_gff        : genome_info[5],
            full_length_TE_gff: genome_info[8],
            full_length_copies: genome_info[9],
            intact_LTR_list: genome_info[10].toString()
        ]
    }.collect().map { data ->
        def jsonContent = "[\n" + data.collect { JsonOutput.toJson(it) }.join(",\n") + "\n]"
        def filePath = "${params.out_dir}/genome_info.json"
        new File(filePath).write(jsonContent)
        return filePath
    }.set { genome_info_json }

    if (!params.skip_analyze) {
         // Step 6: 对检测到的 TE 进行统计分析
        pan_summarize_tes(genome_info_json, params.pan_genomes_dir, panTE_merge_lib, params.softcore_threshold)

        // Step 7: 基因和 TE 关系
        gene_te_associations_out = pan_gene_te_relation(genome_info_json)

        //Step 8: 为RNA_seq生成比对bam
        genome_info_list.map { genome_name, raw_name, reference, gene_gtf, RNA_seq ->
            [genome_name, reference, params.RNA_dir, RNA_seq]
        }.set { generate_bam_input_channel }
        bam_out = pan_generate_bam_for_RNA_seq(generate_bam_input_channel)

        // Step 9: 检测差异表达基因
        // 将 bam 结果合并到 json 文件中, 为pan_detect_de_genes生成输入channel
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
                genome_info[7].toString(),
                genome_info[8].toString(),
                genome_info[9].toString()
            ]
        }.join(bam_out).map { genome_info ->
            [
                genome_name: genome_info[0],
                raw_name: genome_info[1],
                reference: genome_info[2],
                gene_gtf: genome_info[3],
                RNA_seq: genome_info[4],
                TE_gff: genome_info[5],
                full_length_TE_gff: genome_info[8],
                bam: genome_info[10].toString()
            ]
        }.collect().map { data ->
            def jsonContent = "[\n" + data.collect { JsonOutput.toJson(it) }.join(",\n") + "\n]"
            def filePath = "${params.out_dir}/genome_info_for_bam.json"
            new File(filePath).write(jsonContent)
            return filePath
        }.set { genome_info_for_bam_json }
        // genome_name, raw_name, reference, gene_gtf, RNA_seq, TE_gff, TE_tbl, TE_out, TE_full_length_gff, TE_full_length_copies, bam

        pan_detect_de_genes(genome_info_for_bam_json, gene_te_associations_out, params.RNA_dir)
    }
}