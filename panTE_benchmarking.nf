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

// 定义用户输入参数
params.pan_genomes_dir = ''
params.genome_list = ''
params.TE_lib = ''
params.EDTA_home = ''
params.species = ''
params.out_dir = './output'
params.threads = 10

// Step 3: HiTE 并行处理每个基因组
process run_benchmarking_single {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/run_benchmarking_single/${genome_name}"

    input:
    val genome_name

    output:
    path "BM_HiTE.log",    emit: ch_BM_HiTE
    path "BM_EDTA.log",    emit: ch_BM_EDTA

    script:
    cores = task.cpus
    """
    benchmarking.py --BM_EDTA 1 --BM_HiTE 1 -t ${cores} --TE_lib ${params.TE_lib} \
    -r ${params.pan_genomes_dir}/${genome_name} --species ${params.species} --EDTA_home ${params.EDTA_home} --recover 1
    """
}



// 定义工作流
workflow {
    // 读取文件并创建一个 Channel
    Channel.fromPath(params.genome_list, type: 'any', checkIfExists: true)
    .splitText()  // 将文件按行分割
    .map { it.trim() }  // 去掉每行前后的空白字符和换行符
    .filter { it != "" }  // 过滤掉空行
    .set { genome_names }  // 存储为 genome_names 变量

    // 并行处理每个基因组的评测
    run_benchmarking_single(genome_names)
}
