#!/usr/bin/env python
import argparse
import csv
import os
import sys
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import rename_reference, file_exist, Logger, run_LTR_detection, run_LTR_retriever, \
    read_fasta, store_fasta, run_HybridLTR, assign_label_to_lib


def rename_LTR(ltr_lib, rename_ltr_lib):
    contigNames, contigs = read_fasta(ltr_lib)
    chr_elements = {}
    for i, name in enumerate(contigNames):
        raw_name = name.split('#')[0]
        parts = raw_name.split(':')
        chr_name = parts[0]
        info_parts = parts[1].split('_')
        pos_parts = info_parts[0].split('..')
        ltr_type = info_parts[1]
        start_pos = int(pos_parts[0])
        end_pos = int(pos_parts[1])
        element = (start_pos, end_pos, ltr_type, name)
        if not chr_elements.__contains__(chr_name):
            chr_elements[chr_name] = []
        elements = chr_elements[chr_name]
        elements.append(element)
    ltr_index = 0
    ltrName2ltrIndex = {}
    for chr_name in chr_elements.keys():
        elements = chr_elements[chr_name]
        elements.sort(key=lambda x: (x[0], x[1]))
        for i in range(len(elements)):
            cur_element = elements[i]
            cur_ltr_name = cur_element[3]
            if not ltrName2ltrIndex.__contains__(cur_ltr_name):
                ltrName2ltrIndex[cur_ltr_name] = ltr_index
                ltr_index += 1
            cur_ltr_index = ltrName2ltrIndex[cur_ltr_name]
            for j in range(i + 1, len(elements)):
                next_element = elements[j]
                # 是同一个LTR
                if cur_element[2] != next_element[2] and abs(next_element[0] - cur_element[1]) < 10:
                    next_ltr_name = next_element[3]
                    if not ltrName2ltrIndex.__contains__(next_ltr_name):
                        ltrName2ltrIndex[next_ltr_name] = cur_ltr_index
    rename_LTRs = {}
    for ltrName in ltrName2ltrIndex.keys():
        ltr_type = ltrName.split('#')[0].split('_')[-1]
        ltr_label = ltrName.split('#')[1]
        ltr_index = ltrName2ltrIndex[ltrName]
        new_name = 'ltr_' + str(ltr_index) + '-' + ltr_type + '#' + ltr_label
        rename_LTRs[new_name] = contigs[ltrName]
    store_fasta(rename_LTRs, rename_ltr_lib)

def get_intact_ltr(genome_path, ltr_list, intact_LTR_path):
    segmentLTRs = {}
    ltr_contigs = {}
    ref_contigNames, ref_contigs = read_fasta(genome_path)
    ltr_index = 1
    with open(ltr_list, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            direct = line.split('\t')[8]
            ltr_name = line.split('\t')[0]
            internal_info = line.split('\t')[6]
            int_start = int(internal_info.split(':')[1].split('..')[0])
            int_end = int(internal_info.split(':')[1].split('..')[1])
            parts = ltr_name.split(':')
            chr_name = parts[0]
            chr_start = int(parts[1].split('..')[0])
            chr_end = int(parts[1].split('..')[1])
            if direct == '-':
                temp = chr_end
                chr_end = chr_start
                chr_start = temp
            chr_seq = ref_contigs[chr_name]
            ltr_seq = chr_seq[chr_start - 1: chr_end]

            terminal_seq1 = chr_seq[chr_start - 1: int_start - 1]
            internal_seq = chr_seq[int_start - 1: int_end]
            terminal_seq2 = chr_seq[int_end: chr_end]
            current_name = 'ltr_' + str(ltr_index)
            segmentLTRs[current_name + '_LTR'] = terminal_seq1
            segmentLTRs[current_name + '_INT'] = internal_seq
            ltr_contigs[ltr_name] = ltr_seq
            ltr_index += 1
    store_fasta(ltr_contigs, intact_LTR_path)
    return segmentLTRs

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE LTR module...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='Input genome assembly path')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('--use_HybridLTR', metavar='use_HybridLTR',
                        help='Whether to use HybridLTR to identify LTRs, 1: true, 0: false.')
    parser.add_argument('--use_NeuralTE', metavar='use_NeuralTE',
                        help='Whether to use NeuralTE to classify TEs, 1: true, 0: false.')
    parser.add_argument('--miu', metavar='miu',
                        help='The neutral mutation rate (per bp per ya)')
    parser.add_argument('--is_wicker', metavar='is_wicker',
                        help='Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker.')
    parser.add_argument('--is_output_lib', metavar='is_output_lib',
                        help='Whether to output LTR library.')
    parser.add_argument('--debug', metavar='recover',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false.')


    args = parser.parse_args()
    reference = args.g

    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    recover = args.recover
    use_HybridLTR = int(args.use_HybridLTR)
    use_NeuralTE = int(args.use_NeuralTE)
    miu = args.miu
    is_wicker = args.is_wicker
    is_output_lib = int(args.is_output_lib)
    debug = args.debug

    LTR_harvest_parallel_Home = cur_dir + '/bin/LTR_HARVEST_parallel'
    LTR_finder_parallel_Home = cur_dir + '/bin/LTR_FINDER_parallel-master'
    NeuralTE_home = cur_dir + '/bin/NeuralTE'
    TEClass_home = cur_dir + '/classification'
    HybridLTR_home = cur_dir + '/bin/HybridLTR-main'

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir + '/HiTE_ltr.log', level='debug')

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    if debug is None:
        debug = 0
    else:
        debug = int(debug)

    # rename reference
    ref_rename_path = tmp_output_dir + '/genome.rename.fa'
    chr_name_map = tmp_output_dir + '/chr_name.map'
    rename_reference(reference, ref_rename_path, chr_name_map)

    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    confident_ltr_terminal = tmp_output_dir + '/confident_ltr.terminal.fa'
    confident_ltr_internal = tmp_output_dir + '/confident_ltr.internal.fa'
    confident_intact_ltr_list = tmp_output_dir + '/intact_LTR.list'
    intact_LTR_path = tmp_output_dir + '/intact_LTR.fa'

    if use_HybridLTR:
        LTR_output_dir = tmp_output_dir + '/HybridLTR_output'
        confident_ltr = LTR_output_dir + '/confident_ltr.fa'
        Hybrid_ltr_terminal = LTR_output_dir + '/confident_ltr.terminal.fa'
        Hybrid_ltr_internal = LTR_output_dir + '/confident_ltr.internal.fa'
        Hybrid_intact_ltr_list = LTR_output_dir + '/intact_LTR.list'
        HybridLTR_intact_LTR_path = LTR_output_dir + '/intact_LTR.fa'

        check_files = [
            confident_ltr_terminal, confident_ltr_internal, confident_ltr_cut_path,
            confident_intact_ltr_list, intact_LTR_path
        ]

        # 判断是否需要重新运行 LTR 模块
        is_rerun = not all(file_exist(f) for f in check_files)
        if not is_recover or is_rerun:
            run_HybridLTR(ref_rename_path, LTR_output_dir, HybridLTR_home, threads, miu, recover, debug, is_output_lib,
                          log)
            if file_exist(HybridLTR_intact_LTR_path):
                os.system('cp ' + HybridLTR_intact_LTR_path + ' ' + intact_LTR_path)
            if file_exist(Hybrid_ltr_terminal):
                os.system('cp ' + Hybrid_ltr_terminal + ' ' + confident_ltr_terminal)
            if file_exist(Hybrid_ltr_internal):
                os.system('cp ' + Hybrid_ltr_internal + ' ' + confident_ltr_internal)
            if file_exist(Hybrid_intact_ltr_list):
                os.system('cp ' + Hybrid_intact_ltr_list + ' ' + confident_intact_ltr_list)
            if file_exist(confident_ltr):
                os.system('cp ' + confident_ltr + ' ' + confident_ltr_cut_path)
            else:
                log.logger.info(
                    'No LTR retrotransposons are detected in the genome, HiTE continues to identify other types of transposons')
        else:
            for check_file in check_files:
                log.logger.info(check_file + ' exists, skip...')

    else:
        check_files = [
            confident_ltr_cut_path, confident_intact_ltr_list, intact_LTR_path
        ]
        # 判断是否需要重新运行 LTR 模块
        is_rerun = not all(file_exist(f) for f in check_files)
        if not is_recover or is_rerun:
            ltrharvest_output = ref_rename_path + '.harvest.combine.scn'
            ltrfinder_output = ref_rename_path + '.finder.combine.scn'
            if not is_recover or not file_exist(ltrharvest_output) or not file_exist(ltrfinder_output):
                starttime = time.time()
                log.logger.info('Start step0.1: Running LTR_harvest_parallel and LTR_finder_parallel')
                run_LTR_detection(ref_rename_path, tmp_output_dir, threads, LTR_harvest_parallel_Home,
                                  LTR_finder_parallel_Home, log)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step0.1: %.8s s" % (dtime))
            else:
                log.logger.info(ltrharvest_output + ' exists, skip...')
                log.logger.info(ltrfinder_output + ' exists, skip...')

            ltrharvest_output = ref_rename_path + '.harvest.combine.scn'
            ltrfinder_output = ref_rename_path + '.finder.combine.scn'
            ltr_output = tmp_output_dir + '/genome_all.fa.rawLTR.scn'
            os.system('cat ' + ltrharvest_output + ' ' + ltrfinder_output + ' > ' + ltr_output)

            starttime = time.time()
            log.logger.info('Start step0.2: run LTR_retriever to get confident LTR')
            confident_ltr = ref_rename_path + '.LTRlib.fa'
            run_LTR_retriever(ref_rename_path, tmp_output_dir, threads, miu, log)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step0.2: %.8s s" % (dtime))

            # get intact LTR from LTR_retriever pass list
            ltr_list = ref_rename_path + '.pass.list'
            if os.path.exists(ltr_list):
                segmentLTRs = get_intact_ltr(ref_rename_path, ltr_list, intact_LTR_path)

                # 1. recover intact-LTRs from 'ltr_cons_path'
                # 2. classify intact-LTRs and assign label to 'ltr_cons_path'
                ltr_names, ltr_contigs = read_fasta(confident_ltr_cut_path)
                intact_ltr_names, intact_ltr_contigs = read_fasta(intact_LTR_path)
                filter_intact_ltr_contigs = {}
                for name in ltr_names:
                    seq = ltr_contigs[name]
                    name = name.split('#')[0]
                    intact_ltr_name = name[:-4]
                    if intact_ltr_contigs.__contains__(intact_ltr_name):
                        filter_intact_ltr_contigs[intact_ltr_name] = intact_ltr_contigs[intact_ltr_name]
                    else:
                        filter_intact_ltr_contigs[intact_ltr_name] = seq
                store_fasta(filter_intact_ltr_contigs, intact_LTR_path)
            if file_exist(ltr_list):
                os.system('cp ' + ltr_list + ' ' + confident_intact_ltr_list)
            if file_exist(confident_ltr):
                os.system('cp ' + confident_ltr + ' ' + confident_ltr_cut_path)
            else:
                log.logger.info('No LTR retrotransposons are detected in the genome, HiTE continues to identify other types of transposons')
        else:
            for check_file in check_files:
                log.logger.info(check_file + ' exists, skip...')


    # 初始化文件检查列表
    classified_TE_path = intact_LTR_path + '.classified'
    # classify intact-LTRs
    if use_NeuralTE:
        # classify LTR using NeuralTE
        NeuralTE_output_dir = tmp_output_dir + '/NeuralTE_LTR'
        if not os.path.exists(NeuralTE_output_dir):
            os.makedirs(NeuralTE_output_dir)
        NeuralTE_command = 'python ' + NeuralTE_home + '/src/Classifier.py --data ' + intact_LTR_path \
                           + ' --use_TSD 0 --model_path ' \
                           + NeuralTE_home + '/models/NeuralTE_model.h5 --out_dir ' \
                           + NeuralTE_output_dir + ' --thread ' + str(threads) + ' --is_wicker ' + str(is_wicker)
        log.logger.debug(NeuralTE_command)
        os.system(NeuralTE_command)
        NeuralTE_output = NeuralTE_output_dir + '/classified_TE.fa'
        os.system('cp ' + NeuralTE_output + ' ' + classified_TE_path)
    else:
        # classify LTR using RepeatClassifier
        sample_name = 'test'
        TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                          + ' --consensus ' + intact_LTR_path + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(
            threads) + ' -o ' + tmp_output_dir
        log.logger.debug(TEClass_command)
        os.system(TEClass_command)


    # assign intact LTR labels to `genome.rename.fa.LTRlib.fa`
    classified_names, classified_contigs = read_fasta(classified_TE_path)
    intact_LTR_labels = {}
    for name in classified_names:
        parts = name.split('#')
        intact_LTR_labels[parts[0]] = parts[1]

    assign_label_to_lib(confident_ltr_cut_path, intact_LTR_labels)
    assign_label_to_lib(confident_ltr_terminal, intact_LTR_labels)
    assign_label_to_lib(confident_ltr_internal, intact_LTR_labels)

    # add classified labels to intact_ltr.list
    if os.path.exists(confident_intact_ltr_list):
        temp_file = confident_intact_ltr_list + '.temp.tsv'
        with open(confident_intact_ltr_list, 'r') as infile, open(temp_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')

            for row in reader:
                if not row[0].startswith('#'):
                    # 将分类标签添加到 LTR_loc 对应的行
                    LTR_loc = row[0]
                    label = intact_LTR_labels.get(LTR_loc, 'NA')  # 如果没有找到分类标签，则写入 'NA'

                    # 写入新的行：在倒数第二列插入标签列
                    writer.writerow(row[:-1] + [label] + row[-1:])
        # 替换原始文件为带有分类标签的新文件
        os.replace(temp_file, confident_intact_ltr_list)




