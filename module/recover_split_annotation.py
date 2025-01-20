#!/usr/bin/env python
import sys

def recover_split_annotation(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                # outfile.write(line)  # 保留注释行
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                outfile.write(line)  # 保留格式不正确的行
                continue

            # 处理第一列
            chrom_index = parts[0].split("$")
            if len(chrom_index) == 1:
                # 如果没有 $，则直接使用 chrom
                chrom = chrom_index[0]
                index = 0
            else:
                # 取 $ 前面的部分为 chrom，后面的部分为 index
                chrom = "$".join(chrom_index[:-1])  # 合并 $ 前面的部分
                index = int(chrom_index[-1])  # 取最后一个部分为 index

            # 更新第 4 列和第 5 列
            parts[0] = chrom
            parts[3] = str(int(parts[3]) + index)
            parts[4] = str(int(parts[4]) + index)

            # 写回文件
            outfile.write("\t".join(parts) + "\n")

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    recover_split_annotation(input_file, output_file)