#!/usr/bin/env python
import argparse
import os
import re
import shutil
import sys
import uuid

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import Logger, create_or_clear_directory, copy_files, clean_old_tmp_files_by_dir


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
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    args = parser.parse_args()
    threads = int(args.t)
    classified_TE_consensus = args.classified_TE_consensus
    tmp_output_dir = args.tmp_output_dir
    annotate = args.annotate
    reference = args.r
    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)

    classified_TE_consensus = os.path.abspath(classified_TE_consensus)
    reference = os.path.abspath(reference)

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()
    tmp_output_dir = os.path.abspath(tmp_output_dir)

    log = Logger(tmp_output_dir+'/HiTE_annotate_genome.log', level='debug')

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'annotate_genome_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)

        annotate_genome(temp_dir, classified_TE_consensus, reference, annotate, threads, log)

        # 计算完之后将结果拷贝回输出目录
        copy_files(temp_dir, tmp_output_dir)

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