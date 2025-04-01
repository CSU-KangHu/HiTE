#!/usr/bin/env python
import argparse
import shutil
import uuid
import os
import subprocess
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)

from Util import clean_old_tmp_files_by_dir, create_or_clear_directory, read_fasta


def run_minimap2(assembly_file, short_fasta, output_paf, threads):
    """Run minimap2 for self-alignment with full contig names."""
    # Run minimap2 with specified threads and output to the temporary directory
    command = f"minimap2 -x asm5 -t {threads} {assembly_file} {short_fasta} > {output_paf}"
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True, check=True)


def filter_redundant_contigs(paf_file, output_txt):
    """Filter redundant contigs using awk and sort."""
    command = f"awk '$1 != $6 && $10/$2 > 0.95 {{print $1}}' {paf_file} | sort | uniq > {output_txt}"
    subprocess.run(command, shell=True, check=True)


def remove_redundant_contigs(assembly_file, redundant_txt, output_assembly):
    """Remove redundant contigs from the assembly."""
    with open(redundant_txt, 'r') as f:
        redundant_contigs = set(f.read().splitlines())

    with open(assembly_file, 'r') as infile, open(output_assembly, 'w') as outfile:
        write = False
        for line in infile:
            if line.startswith('>'):
                contig_id = line.strip().split()[0][1:]  # Extract contig ID
                write = contig_id not in redundant_contigs
            if write:
                outfile.write(line)

def extract_small_sequences(input_fasta, output_fasta, max_size=10000000):
    """
    从输入的FASTA文件中提取小于指定大小的序列，并保存到新的FASTA文件中。
    :param input_fasta: 输入的FASTA文件路径
    :param output_fasta: 输出的FASTA文件路径
    :param max_size: 最大序列大小（默认10MB = 10,000,000碱基对）
    """
    contig_names, contigs = read_fasta(input_fasta)
    with open(output_fasta, "w") as output_handle:
        for name in contig_names:
            seq = contigs[name]
            if len(seq) < max_size:
                output_handle.write(f">{name}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Remove redundant contigs from a genome assembly using minimap2.")
    parser.add_argument('-i', '--input', required=True, help="Input genome assembly file (FASTA format).")
    parser.add_argument('-o', '--output', required=True, help="Output genome assembly file (FASTA format).")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads to use for minimap2.")
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    args = parser.parse_args()
    input_path = args.input
    input_path = os.path.abspath(input_path)
    output_path = args.output
    output_path = os.path.abspath(output_path)

    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'genome_clean_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)

        # Define paths for intermediate files
        paf_file = os.path.join(temp_dir, 'genome_self.asm5.paf')
        redundant_txt = os.path.join(temp_dir, 'redundant_contigs.txt')
        short_fasta = os.path.join(temp_dir, 'short.fasta')

        extract_small_sequences(input_path, short_fasta, max_size=10000000)

        # Step 1: Run minimap2 for self-alignment
        run_minimap2(input_path, short_fasta, paf_file, args.threads)

        # Step 2: Filter redundant contigs
        filter_redundant_contigs(paf_file, redundant_txt)

        # Step 3: Remove redundant contigs and generate new assembly
        remove_redundant_contigs(input_path, redundant_txt, output_path)
    except Exception as e:
        # 如果出现异常，打印错误信息并删除临时目录
        print(f"An error occurred: {e}")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        raise  # 重新抛出异常，以便上层代码可以处理

    else:
        # 如果没有异常，删除临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

if __name__ == "__main__":
    main()




