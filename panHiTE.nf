import groovy.json.JsonSlurper
import groovy.json.JsonOutput

// 定义用户输入参数
params.pan_genomes_dir = ''
params.genome_list = ''
params.genes_dir = ''
params.RNA_dir = ''
params.out_dir = './output'
params.te_type = 'ltr'
params.skip_analyze = false
params.recover = false
params.debug = false
params.threads = 4
params.miu = 0.1
params.all_te_types = ['ltr', 'tir', 'helitron', 'non-ltr', 'all']

// 验证 TE 类型是否合法
if (!params.all_te_types.contains(params.te_type)) {
    error "Invalid TE type: ${params.te_type}. Please choose from ${params.all_te_types}"
}


process preprocess_genomes {
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
    input:
    tuple val(genome_name), val(raw_name), path(reference), val(threads), val(te_type), val(miu), val(debug), val(recover)

    output:
    tuple val(genome_name), path("intact_LTR.list"), emit: ch_intact_ltr_list, optional: true
    path "confident_ltr.terminal.fa", emit: ch_ltr_terminal, optional: true
    path "confident_ltr.internal.fa", emit: ch_ltr_internal, optional: true
    path "confident_helitron.fa", emit: ch_helitron, optional: true
    path "confident_non_ltr.fa", emit: ch_non_ltr, optional: true
    path "confident_other.fa", emit: ch_other, optional: true
    path "confident_tir.fa", emit: ch_tir, optional: true
    path "confident_TE.cons.fa", emit: ch_te, optional: true

    script:
    """
    pan_run_hite_single.py --genome_name ${genome_name} --reference ${reference} --threads ${threads} \
    --te_type ${te_type} --miu ${miu} --debug ${debug} --recover ${recover}
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
    input:
    tuple val(genome_name), path(reference), val(output_dir), val(threads), val(recover), path(panTE_lib)

    output:
    tuple val(genome_name), path("${output_dir}/${genome_name}.gff"), path("${output_dir}/${genome_name}.full_length.gff"), path("${output_dir}/${genome_name}.full_length.copies"), emit: annotate_out

    script:
    """
    pan_annotate_genome.py --threads ${threads} --panTE_lib ${panTE_lib} --reference ${reference} \
    --genome_name ${genome_name} --recover ${recover} --output_dir ${output_dir}
    """
}

// Step 4: 汇总 TE 数据
process summarize_tes {
    input:
    path genome_info_json
    path pan_genomes_dir
    path panTE_lib
    val recover

    output:
    path "TE_summary.pdf", emit: ch_te_summary
    path "panHiTE.CorePan_fitmodel.pdf", emit: ch_corepan_model, optional: true
    path "panHiTE.CorePan_fitsmooth.pdf", emit: ch_corepan_smodel, optional: true

    script:
    """
    pan_summary_TEs.py --genome_info_json ${genome_info_json} --pan_genomes_dir ${pan_genomes_dir} \
    --panTE_lib ${panTE_lib} --recover ${recover}
    """
}


process pan_gene_te_relation {
    input:
    path genome_info_json
    val recover

    output:
    path "gene_te_associations.tsv"

    script:
    """
    pan_gene_te_relation.py --genome_info_json ${genome_info_json} --recover ${recover}
    """
}

process pan_detect_de_genes {
    input:
    path genome_info_json
    path gene_te_associations
    val RNA_dir
    val threads
    val recover

    output:
    path "DE_genes_from_TEs.tsv", emit: ch_de_genes, optional: true
    path "all_gene_TEs_details.tsv", emit: ch_all_genes, optional: true

    script:
    """
    pan_detect_de_genes.py --genome_info_json ${genome_info_json} --gene_te_associations ${gene_te_associations} \
    --RNA_dir ${RNA_dir} --threads ${threads} --recover ${recover}
    """
}


// 定义工作流
workflow {
    // Step 1: 预处理，生成json格式的输入文件路径
    genome_metadata_out = preprocess_genomes(params.genome_list, params.genes_dir, params.RNA_dir, params.pan_genomes_dir)
    genome_metadata_out.collectFile(name: "${params.out_dir}/genome_metadata.json")

    // Step 2: 解析Step1的json文件，准备HiTE输入channel
    genome_info_list = genome_metadata_out
        .flatMap { json_file ->
            def genome_data = file(json_file).text  // 读取文件内容为字符串
            def genome_json = new JsonSlurper().parseText(genome_data) // 使用 JsonSlurper 解析 JSON
            genome_json.genome_info.collect { [it.genome_name, it.raw_name, it.reference, it.gene_gtf, it.RNA_seq] }
        }
    genome_info_list.map { genome_name, raw_name, reference, gene_gtf, RNA_seq ->
            [genome_name, raw_name, reference, params.threads, params.te_type, params.miu, params.debug, params.recover]
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
        [genome_name, reference, params.out_dir, params.threads, params.recover]
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
        te_summary_out = summarize_tes(genome_info_json, params.pan_genomes_dir, panTE_lib, params.recover)
        te_summary_out.ch_te_summary.collectFile(name: "${params.out_dir}/TE_summary.pdf")
        te_summary_out.ch_corepan_model.collectFile(name: "${params.out_dir}/panHiTE.CorePan_fitmodel.pdf")
        te_summary_out.ch_corepan_smodel.collectFile(name: "${params.out_dir}/panHiTE.CorePan_fitsmooth.pdf")

        // Step 8: 基因和 TE 关系
        gene_te_associations_out = pan_gene_te_relation(genome_info_json, params.recover)
        gene_te_associations_out.collectFile(name: "${params.out_dir}/gene_te_associations.tsv")

        // Step 9: 检测差异表达基因
        detected_de_genes_out = pan_detect_de_genes(genome_info_json, gene_te_associations_out, params.RNA_dir, params.threads, params.recover)
        detected_de_genes_out.ch_de_genes.collectFile(name: "${params.out_dir}/DE_genes_from_TEs.tsv")
        detected_de_genes_out.ch_all_genes.collectFile(name: "${params.out_dir}/all_gene_TEs_details.tsv")
    }
}
