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
    path out_dir

    output:
    path "${out_dir}/genome_metadata.json"

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
    path "HiTE_${raw_name}/${genome_name}_hite_result.json"

    script:
    """
    pan_run_hite_single.py --genome_name ${genome_name} --reference ${reference} --threads ${threads} \
    --te_type ${te_type} --miu ${miu} --debug ${debug} --recover ${recover}
    """
}

// Step 4: 合并 HiTE 结果
process merge_hite_results {
    input:
    file hite_results  // 传递 HiTE 结果文件作为输入

    output:
    path "pan_terminal.tmp.fa",    emit: ch_terminal, optional: true
    path "pan_internal.tmp.fa",    emit: ch_internal, optional: true
    path "intact_ltr_paths.txt",    emit: ch_intact_ltr, optional: true

    script:
    """
    // 创建临时文件，存储合并的 LTR 信息
    terminal_tmp=\${output_dir}/pan_terminal.tmp.fa
    internal_tmp=\${output_dir}/pan_internal.tmp.fa

    // 清空临时文件
    echo -n > \${terminal_tmp}
    echo -n > \${internal_tmp}

    // 创建一个 list 用于保存所有基因组的 intact LTR 路径
    intact_ltr_paths=()

    // 遍历 HiTE 结果文件并提取需要的字段
    for result_file in ${hite_results}; do
        // 读取当前结果文件
        single_result=\$(cat \${result_file})

        // 提取所需的文件路径
        confident_ltr_terminal=\$(echo \${single_result} | jq -r '.confident_ltr_terminal')
        confident_ltr_internal=\$(echo \${single_result} | jq -r '.confident_ltr_internal')
        confident_TE=\$(echo \${single_result} | jq -r '.confident_TE')
        ltr_intact_list=\$(echo \${single_result} | jq -r '.ltr_intact_list')

        // 将 LTR 路径添加到 intact_ltr_paths 列表
        intact_ltr_paths+=(\${ltr_intact_list})

        // 将各个结果文件的内容合并到临时文件
        cat \${confident_ltr_terminal} >> \${terminal_tmp}
        cat \${confident_ltr_internal} >> \${internal_tmp}
        cat \${confident_TE} >> \${terminal_tmp}
    done

    // 生成 intact_ltr_paths_file 并保存路径信息
    echo \$(echo \${intact_ltr_paths} | jq -c .) > \${output_dir}/intact_ltr_paths.txt
    """
}

// Step 5: 去冗余
process pan_remove_redundancy {
    input:
    file terminal_tmp
    file internal_tmp
    file intact_ltr_paths_file
    val out_dir
    val threads

    output:
    file("${out_dir}/panTE.fa")

    script:
    """
    if (params.recover && file("${out_dir}/panTE.fa").exists()) {
        println "Skipping redundancy removal as panTE.fa already exists and recover is enabled"
    } else {
        println "Running redundancy removal"
        pan_remove_redundancy.py ${terminal_tmp} ${internal_tmp} ${out_dir} ${threads} ${out_dir}/panTE.fa
    }
    """
}

// Step 3: 注释基因组
process annotate_genomes {
    input:
    tuple val(genome_name), path(reference), path(TE_gff), path(panTE_lib), val(out_dir), val(threads), val(recover)


    output:
    file("${out_dir}/${genome_name}_annotation.json")

    script:
    """
    pan_annotate_genome.py ${TE_gff} ${out_dir} ${threads} ${panTE_lib} ${reference} ${genome_name} ${recover}
    """
}

// Step 4: 汇总 TE 数据
process summarize_tes {
    input:
    file genome_metadata
    file pan_genomes_dir
    file panTE_lib
    file intact_ltr_paths_file
    file total_annotation
    val recover

    output:
    file("${params.out_dir}/summary_TEs.json")

    script:
    """
    pan_summary_TEs.py ${genome_metadata} ${pan_genomes_dir} ${panTE_lib} ${intact_ltr_paths_file} ${recover}
    """
}


process pan_gene_te_relation {
    input:
    file genome_metadata
    file total_annotation
    val out_dir, recover

    output:
    file("${out_dir}/gene_te_associations.tsv")

    script:
    """
    if (params.recover && file("${out_dir}/gene_te_associations.tsv").exists()) {
        echo "Skipping gene-te relation as the result already exists and recover is enabled"
    } else {
        echo "Running gene-te relation analysis"
        gene_te_associations="${out_dir}/gene_te_associations.tsv"
        pan_gene_te_relation.py ${genome_metadata} ${out_dir} ${recover}
    }
    """
}

process pan_detect_de_genes {
    input:
    file genome_metadata
    file gene_te_associations
    val out_dir, threads, recover

    output:
    file("${out_dir}/detected_de_genes.json")

    script:
    """
    if (params.recover && file("${out_dir}/detected_de_genes.json").exists()) {
        echo "Skipping DE gene detection as the result already exists and recover is enabled"
    } else {
        echo "Running DE gene detection"
        bash pan_detect_de_genes.py ${genome_metadata} ${threads} ${recover} ${out_dir} ${gene_te_associations}"
    }
    """
}

import groovy.json.JsonSlurper

// 定义工作流
workflow {
    // Step 1: 预处理
    //genome_metadata_out = preprocess_genomes(params.genome_list, params.genes_dir, params.RNA_dir, params.pan_genomes_dir, params.out_dir)

    Channel.fromPath('/home/hukang/test/HiTE/demo/demo/panHiTE_output_demo_nf/genome_metadata.json', type: 'any', checkIfExists: true).set{ genome_metadata_out }

    // Step 2: 从 genome_metadata.json 加载 genome_info_list
    genome_info_list = genome_metadata_out
        .flatMap { json_file ->
            def genome_data = file(json_file).text  // 读取文件内容为字符串
            def genome_json = new JsonSlurper().parseText(genome_data) // 使用 JsonSlurper 解析 JSON
            genome_json.genome_info.collect { [it.genome_name, it.raw_name, it.reference] }
        }

    // Step 3: HiTE 并行处理每个基因组
    hite_input = genome_info_list
        .map { genome_name, raw_name, reference ->
            [genome_name, raw_name, reference, params.threads, params.te_type, params.miu, params.debug, params.recover]
        }
        .set { hite_input_channel }

    hite_out = run_hite_single(hite_input_channel)
//
//     // Step 4: 合并 HiTE 结果（提取文件并合并）
//     merged_hite_files = merge_hite_results(hite_out)
//
//     // 获取 intact_ltr_paths_file
//     def intact_ltr_paths_file = merged_hite_files.ch_intact_ltr
//
//     // Step 5: 去冗余
//     panTE_lib = pan_remove_redundancy(merged_hite_files, params.out_dir, params.threads)
//
//     // Step 6: 注释基因组
//     annotate_input = genome_info_list
//         .map { genome_name, reference, TE_gff ->
//             [genome_name, reference, TE_gff, panTE_lib, params.out_dir, params.threads, params.recover]
//         }
//         .set { annotate_input_channel }
//
//     annotate_out = annotate_genomes(annotate_input_channel)
//
//     total_annotation = annotate_out.collectFile(name: "${params.out_dir}/total_annotation.json")
//
//     if (!params.skip_analyze) {
//          // Step 7: 汇总 TE 数据
//         summarize_tes(genome_metadata_out, params.pan_genomes_dir, panTE_lib, intact_ltr_paths_file, total_annotation, params.recover)
//
//         // Step 8: 基因和 TE 关系
//         gene_te_associations_out = pan_gene_te_relation(genome_metadata_out, total_annotation, params.out_dir, params.recover)
//
//         // Step 9: 检测差异表达基因
//         detected_de_genes_out = pan_detect_de_genes(genome_metadata_out, gene_te_associations_out, params.out_dir, params.threads, params.recover)
//     }
}
