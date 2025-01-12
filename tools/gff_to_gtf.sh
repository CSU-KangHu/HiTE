#!/bin/bash

# 输入 GFF3 文件
input_gff3=$1

# 输出 GTF 文件
output_gtf="${input_gff3%.*}.gtf"

# 转换函数
awk '
BEGIN {
    FS="\t";  # 输入字段分隔符
    OFS="\t"; # 输出字段分隔符
}

# 跳过注释行
/^#/ { next; }

{
    # 提取前 8 个字段
    seqid = $1;
    source = $2;
    type = "exon";
    start = $4;
    end = $5;
    score = $6;
    strand = $7;
    phase = $8;

    # 初始化 gene_id 和 transcript_id
    gene_id = "unknown";
    transcript_id = "unknown";

    # 处理 Target 属性格式
    if (match($9, /Target "([^"]+)"/, target_match)) {
        target_value = target_match[1];  # 提取 "Motif:19-LTR_253-int"
        gsub(/Motif:/, "", target_value);  # 移除 "Motif:"
        gene_id = target_value;
        transcript_id = target_value;
    }
    # 处理 id= 和 name= 属性格式
    else if (match($9, /id=[^;]+/)) {
        id_value = substr($9, RSTART + 3, RLENGTH - 3);  # 提取 "te_intact_4"
        if (match($9, /name=[^;]+/)) {
            name_value = substr($9, RSTART + 5, RLENGTH - 5);  # 提取 "44-Helitron_0"
            gene_id = name_value;
            transcript_id = name_value;
        } else {
            gene_id = id_value;
            transcript_id = id_value;
        }
    }

    # 构造 GTF 属性字段
    attributes = "gene_id \"" gene_id "\";"

    # 打印 GTF 格式
    print seqid, source, type, start, end, score, strand, phase, attributes;
}
' "$input_gff3" > "$output_gtf"

