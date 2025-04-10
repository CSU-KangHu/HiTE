#!/usr/bin/env python
import sys
from datetime import datetime

def load_te_type_mapping(te_type_file_path):
    """加载TE ID和类型的映射关系字典"""
    te_type_mapping = {}
    with open(te_type_file_path, "r") as te_type_file:
        for line in te_type_file:
            if line.strip():
                te_id, te_type = line.strip().split("\t")
                te_type_mapping[te_id] = te_type
    return te_type_mapping

def process_gff(input_file, output_file, te_type_mapping):
    """处理GFF文件，替换第三列的TE ID为TE类型，并输出到输出文件"""
    with open(input_file, "r") as gff_file, open(output_file, "w") as output_file:
        # 获取当前时间
        current_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
        # 写入头部信息
        output_file.write("##gff-version 3\n")
        output_file.write(f"##date {current_time}\n")  # 使用当前时间
        output_file.write("##seqid\tsource\tclassification\tstart\tend\tdivergence\tstrand\tphase\tattributes\n")

        for line in gff_file:
            if line.strip():
                columns = line.strip().split("\t")
                if len(columns) >= 9:
                    # 提取TE ID
                    te_id = columns[8].split("Motif:")[1].split('"')[0]
                    te_type = te_type_mapping.get(te_id, "SSR/low_complexity")
                    columns[1] = 'HiTE'
                    # 替换第三列信息为TE类型
                    columns[2] = te_type
                    # 构建新的attributes
                    attributes = f'Name={te_id};Classification={te_type}'
                    columns[8] = attributes
                    # 写入处理后的行
                    output_file.write("\t".join(columns) + "\n")
                else:
                    # 保留格式不正确的行
                    output_file.write(line)

def main():
    te_type_file_path = sys.argv[1]
    gff_input_file_path = sys.argv[2]
    gff_output_file_path = sys.argv[3]

    # 加载TE类型映射
    te_type_mapping = load_te_type_mapping(te_type_file_path)

    # 处理GFF文件
    process_gff(gff_input_file_path, gff_output_file_path, te_type_mapping)

if __name__ == "__main__":
    main()