#!/usr/bin/env python
import argparse
import gzip
import os
import sys
import json

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import convertGeneAnnotation2GTF, read_fasta, store_fasta, Logger

def preprocess_genomes(genome_list_path, genes_dir, RNA_dir, pan_genomes_dir, output_dir, log):
    genome_paths = []
    genome_names = []
    is_RNA_analyze = False
    gene_annotation_list = []

    # 1.1. 读取 genomes, gff, RNA-seq reads
    with open(genome_list_path, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            genome_name = parts[0]
            gene_gff = None

            if len(parts) >= 2:
                gene_name = parts[1]
                if gene_name:
                    gene_gff = os.path.join(genes_dir, gene_name)
                    if not os.path.isabs(gene_gff):
                        gene_gff = os.path.abspath(gene_gff)
                    if not os.path.exists(gene_gff):
                        log.logger.error(f'Cannot find gene annotation path: {gene_gff}')
                        sys.exit(-1)
                    if not (gene_gff.endswith('.gff') or gene_gff.endswith('.gff3')):
                        log.logger.error(f'Gene annotation file must ends with .gff or .gff3')
                        sys.exit(-1)
                    gene_annotation_list.append(gene_gff)

            RNA_seq_dict = {}
            if len(parts) > 3:
                is_PE = bool(int(parts[2]))
                if is_PE:
                    raw_RNA1 = os.path.join(RNA_dir, parts[3])
                    raw_RNA2 = os.path.join(RNA_dir, parts[4])
                    if not os.path.exists(raw_RNA1) or not os.path.exists(raw_RNA2):
                        log.logger.error(f'Cannot find RNA-seq path: {raw_RNA1}, {raw_RNA2}')
                        sys.exit(-1)
                    RNA_seq_dict['raw_RNA1'] = parts[3]
                    RNA_seq_dict['raw_RNA2'] = parts[4]
                    RNA_seq_dict['is_PE'] = is_PE
                    is_RNA_analyze = True
                else:
                    cur_merge_fq = os.path.join(RNA_dir, genome_name + '.merge.fq.gz')
                    with gzip.open(cur_merge_fq, 'wb') as outfile:
                        for raw_RNA_name in parts[3:]:
                            raw_RNA = os.path.join(RNA_dir, raw_RNA_name)
                            if not os.path.exists(raw_RNA):
                                log.logger.error(f'Cannot find RNA-seq path: {raw_RNA}')
                                sys.exit(-1)
                            # 将文件内容追加到输出文件
                            with gzip.open(raw_RNA, 'rb') as infile:
                                outfile.write(infile.read())
                    RNA_seq_dict['raw_RNA'] = genome_name + '.merge.fq.gz'
                    RNA_seq_dict['is_PE'] = is_PE
                    is_RNA_analyze = True

            cur_genome_path = os.path.join(pan_genomes_dir, genome_name)

            if not os.path.isabs(cur_genome_path):
                cur_genome_path = os.path.abspath(cur_genome_path)
            if not os.path.exists(cur_genome_path):
                log.logger.error(f'Cannot find genome path: {cur_genome_path}')
                sys.exit(-1)

            genome_paths.append((genome_name, cur_genome_path, gene_gff, RNA_seq_dict))
            genome_names.append(genome_name)

    # 1.2. 将 gene 注释文件转为标准的 gtf 格式，并去除非编码基因
    script_dir = os.path.join(project_dir, 'RNA_seq')
    genome_paths = convertGeneAnnotation2GTF(genome_paths, script_dir, output_dir, log)

    # 1.3. 生成 total_genome.fa
    # total_genome = os.path.join(output_dir, 'total_genome.fa')
    # new_ref_contigs = {}
    genome_info_list = []
    for genome_name, reference, gene_gtf, RNA_seq_dict in genome_paths:
        raw_name = genome_name.split('.')[0]
        # ref_names, ref_contigs = read_fasta(reference)
        # for name in ref_names:
        #     new_name = f'{raw_name}-{name}'
        #     new_ref_contigs[new_name] = ref_contigs[name]
        genome_info_list.append({
            "genome_name": genome_name,
            "raw_name": raw_name,
            "reference": reference,
            "gene_gtf": gene_gtf,
            "RNA_seq": RNA_seq_dict
        })
    # store_fasta(new_ref_contigs, total_genome)

    # 保存 JSON 信息
    output_metadata = os.path.join(output_dir, 'genome_metadata.json')
    with open(output_metadata, 'w') as f_out:
        json.dump({
            # "total_genome": total_genome,
            "genome_info": genome_info_list,
            "is_RNA_analyze": is_RNA_analyze
        }, f_out, indent=4)

    log.logger.info('Preprocessing completed.')
    return output_metadata


if __name__ == '__main__':
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE preprocess genomes and gff annotations.")
    parser.add_argument("--genome_list", type=str, help="Path to the genome list file.")
    parser.add_argument("--genes_dir", type=str, help="Path to the directory containing gene files.")
    parser.add_argument("--pan_genomes_dir", type=str, help="Path to the directory containing pan-genomes.")
    parser.add_argument("--RNA_dir", nargs="?", default='', help="Path to the directory containing RNA files.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Path to the output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 接收命令行参数
    genome_list = args.genome_list
    genes_dir = args.genes_dir
    RNA_dir = args.RNA_dir
    pan_genomes_dir = args.pan_genomes_dir

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    preprocess_genomes(genome_list, genes_dir, RNA_dir, pan_genomes_dir, output_dir, log)
