#!/usr/bin/env python
import argparse
import os
import re
import shutil
import sys
import uuid

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import Logger, create_or_clear_directory, copy_files


def annotate_genome(tmp_output_dir, classified_TE_consensus, reference, annotate, threads, log):
    # annotate the genome
    if annotate is not None and int(annotate) == 1:
        RepeatMasker_command = f'cd {tmp_output_dir} && RepeatMasker -e ncbi -pa {threads} -gff -lib {classified_TE_consensus} -cutoff 225 {reference}'
        log.logger.info(f"Running command: {RepeatMasker_command}")
        os.system(RepeatMasker_command)

        mv_file_command = 'mv ' + reference + '.out ' + tmp_output_dir + '/HiTE.out && mv ' \
                          + reference + '.tbl ' + tmp_output_dir + '/HiTE.tbl && mv ' \
                          + reference + '.out.gff ' + tmp_output_dir + '/HiTE.gff && mv ' \
                          + reference + '.cat.gz ' + tmp_output_dir + '/HiTE.cat.gz'
        log.logger.debug(mv_file_command)
        os.system(mv_file_command)

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE annotating genome...')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--classified_TE_consensus', metavar='classified_TE_consensus',
                        help='The path of classified TE library')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--annotate', metavar='annotate',
                        help='Whether to annotate the genome using the TE library generated, 1: true, 0: false.')
    parser.add_argument('-r', metavar='reference',
                        help='Input reference Path')

    args = parser.parse_args()
    threads = int(args.t)
    classified_TE_consensus = args.classified_TE_consensus
    tmp_output_dir = args.tmp_output_dir
    annotate = args.annotate
    reference = args.r

    classified_TE_consensus = os.path.abspath(classified_TE_consensus)
    reference = os.path.abspath(reference)

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()
    tmp_output_dir = os.path.abspath(tmp_output_dir)

    log = Logger(tmp_output_dir+'/HiTE_annotate_genome.log', level='debug')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = '/tmp/annotate_genome_' + str(unique_id)
    create_or_clear_directory(temp_dir)

    annotate_genome(temp_dir, classified_TE_consensus, reference, annotate, threads, log)

    # 计算完之后将结果拷贝回输出目录
    copy_files(temp_dir, tmp_output_dir)

    # 删除临时目录
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
