params.genome_info = "/home/hukang/test/HiTE/demo/genome_list"  // 基因组信息文件
params.chunk_size = 5  // 每个 chunk 的大小
params.chrom_seg_length = 100_000
params.out_dir = "/home/hukang/test/HiTE/demo/results"             // 输出目录
params.threads = 10
params.te_lib = null

// Step 1: 划分基因组
process split_genome {
    storeDir "${params.out_dir}/split_genome/${genome_name}"

    input:
    tuple val(genome_name), path(reference)

    output:
    tuple val(genome_name), path("genome.cut*.fa"), emit: ch_chunks

    script:
    """
    split_genome_chunks.py -g ${reference} --chrom_seg_length ${params.chrom_seg_length} --chunk_size ${params.chunk_size}
    """
}

// Step 2: 并行注释每个 chunk
process annotate_chunk {
    cpus { params.threads ?: 1 }

    storeDir "${params.out_dir}/annotate_chunk/${genome_name}"

    input:
    tuple val(genome_name), path(chunk_fasta), path(panTE_lib)

    output:
    tuple val(genome_name), path("${genome_name}_${ref_index}.gff"), path("${genome_name}_${ref_index}.full_length.gff"), emit: ch_annotated_chunks

    script:
    (full, ref_index) = (chunk_fasta =~ /genome.cut(\d+)\.fa/)[0]
    """
    pan_annotate_genome.py --threads ${task.cpus} --panTE_lib ${panTE_lib} \
    --reference ${chunk_fasta} --genome_name ${genome_name}_${ref_index}
    """
}

// Step 3: 合并每个基因组的注释结果
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
    tuple val(genome_name), path(input_gff), path(full_length_input_gff), path(genome_path), path(te_lib)

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

// 主流程
workflow {
    // 读取基因组信息文件
    genome_info_list = Channel.fromPath(params.genome_info)
        .splitCsv(sep: "\t", header: true)
        .map { row -> [row.genome_name, row.reference] }

    // 准备 panTE 库
    Channel.fromPath(params.te_lib).set{ panTE_merge_lib }

    // 划分基因组
    split_genome_out = split_genome(genome_info_list)

    // 将文件数组展开，并将每个文件与 genome_name 和 panTE_lib 组合
    expanded_chunks = split_genome_out.ch_chunks.flatMap { genome_name, chunk_files ->
        chunk_files.collect { chunk_file ->
            [genome_name, chunk_file]
        }
    }

    expanded_chunks.view { "The value is: $it" }

    // 合并基因组信息和 panTE 库
    annotate_input = expanded_chunks
        .combine(panTE_merge_lib)
        .set { annotate_input_channel }

    // annotate_input_channel.view { "The value is: $it" }

    // 并行注释每个 chunk
    annotate_chunk_out = annotate_chunk(annotate_input_channel)

    // 按基因组名称分组并合并
    annotate_chunk_out.ch_annotated_chunks
        .groupTuple(by: 0)  // 按 genome_name 分组
        .set { grouped_chunks }

    // grouped_chunks.view { "The value is: $it" }

    // 合并每个基因组的注释结果
    merge_annotations_out = merge_annotations(grouped_chunks)
    ch_merge_annotations_out = merge_annotations_out.ch_merged_annotations

    ch_merge_annotations_out.join(genome_info_list).combine(panTE_merge_lib).set { recover_split_annotation_input_channel }

    // 修复 GFF 文件
    recover_split_annotation_out = recover_split_annotation(recover_split_annotation_input_channel)
}
