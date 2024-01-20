#-- coding: UTF-8 --
import argparse
import os
import sys
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import rename_reference, file_exist, Logger, run_LTR_detection, run_LTR_retriever, \
    rename_fasta, deredundant_for_LTR, read_fasta, store_fasta

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
            ltr_contigs[current_name] = ltr_seq
            ltr_index += 1
    store_fasta(ltr_contigs, intact_LTR_path)
    return segmentLTRs

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE LTR module...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='Input genome assembly path')
    parser.add_argument('--ltrharvest_home', metavar='ltrharvest_home',
                        help='Please enter the root directory for ltr_harvest. Use an absolute path.')
    parser.add_argument('--ltrfinder_home', metavar='ltrfinder_home',
                        help='Please enter the root directory for ltr_finder. Use an absolute path.')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('--use_NeuralTE', metavar='use_NeuralTE',
                        help='Whether to use NeuralTE to classify TEs, 1: true, 0: false.')
    parser.add_argument('--NeuralTE_home', metavar='NeuralTE_home',
                        help='The root directory of NeuralTE')
    parser.add_argument('--TEClass_home', metavar='TEClass_home',
                        help='The root directory of TEClass')
    parser.add_argument('--miu', metavar='miu',
                        help='The neutral mutation rate (per bp per ya)')


    args = parser.parse_args()
    reference = args.g
    LTR_harvest_parallel_Home = args.ltrharvest_home
    LTR_finder_parallel_Home = args.ltrfinder_home
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    recover = args.recover
    use_NeuralTE = int(args.use_NeuralTE)
    miu = args.miu
    NeuralTE_home = args.NeuralTE_home
    TEClass_home = args.TEClass_home

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir + '/HiTE_ltr.log', level='debug')

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    # rename reference
    ref_rename_path = tmp_output_dir + '/genome.rename.fa'
    rename_reference(reference, ref_rename_path)

    ltrharvest_output = ref_rename_path + '.harvest.combine.scn'
    ltrfinder_output = ref_rename_path + '.finder.combine.scn'
    if not is_recover or not file_exist(ltrharvest_output) or not file_exist(ltrfinder_output):
        starttime = time.time()
        log.logger.info('Start step0.1: Running LTR_harvest_parallel and LTR_finder_parallel')
        run_LTR_detection(ref_rename_path, tmp_output_dir, threads, LTR_harvest_parallel_Home, LTR_finder_parallel_Home, log)
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

    # resut_file = ref_rename_path + '.LTRlib.fa'
    resut_file = ref_rename_path + '.pass.list'
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        log.logger.info('Start step0.2: run LTR_retriever to get confident LTR')
        run_LTR_retriever(ref_rename_path, tmp_output_dir, threads, miu, log)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step0.2: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')

    # 1. get intact LTR from LTR_retriever pass list
    # 2. get LTR internal and terminal sequences (segment LTR)
    # 3. map segment LTR to intact LTR, so that we can assign labels from intact LTR to segment LTR
    intact_LTR_path = tmp_output_dir + '/intact_LTR.fa'
    segmentLTRs = get_intact_ltr(ref_rename_path, resut_file, intact_LTR_path)

    if use_NeuralTE:
        # classify LTR using NeuralTE
        NeuralTE_output_dir = tmp_output_dir + '/NeuralTE'
        if not os.path.exists(NeuralTE_output_dir):
            os.makedirs(NeuralTE_output_dir)
        NeuralTE_command = 'python ' + NeuralTE_home + '/src/Classifier.py --data ' + intact_LTR_path \
                           + ' --genome ' + ref_rename_path + ' --use_TSD 1 --model_path ' \
                           + NeuralTE_home + '/models/NeuralTE-TSDs_model.h5 --outdir ' \
                           + NeuralTE_output_dir + ' --thread ' + str(threads) + ' --is_wicker 0'
        log.logger.debug(NeuralTE_command)
        os.system(NeuralTE_command + ' > /dev/null 2>&1')
        classified_TE_path = NeuralTE_output_dir + '/classified_TE.fa'
    else:
        # classify LTR using RepeatClassifier
        sample_name = 'test'
        TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                          + ' --consensus ' + intact_LTR_path + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(threads) + ' -o ' + tmp_output_dir
        log.logger.debug(TEClass_command)
        os.system(TEClass_command)
        classified_TE_path = intact_LTR_path + '.classified'


    classified_names, classified_contigs = read_fasta(classified_TE_path)
    intact_LTR_labels = {}
    for name in classified_names:
        parts = name.split('#')
        intact_LTR_labels[parts[0]] = parts[1]

    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    confident_ltr_cut_contigs = {}
    for name in segmentLTRs.keys():
        parts = name.split('_')
        intact_ltr_name = parts[0] + '_' + parts[1]
        label = intact_LTR_labels[intact_ltr_name]
        new_name = name+'#'+label
        confident_ltr_cut_contigs[new_name] = segmentLTRs[name]
    store_fasta(confident_ltr_cut_contigs, confident_ltr_cut_path)

    # rename_LTR(resut_file, confident_ltr_cut_path)
    # confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    # rename_fasta(resut_file, confident_ltr_cut_path, 'LTR')

    # Remove redundancy from the LTR results.
    ltr_cons_path = confident_ltr_cut_path + '.cons'
    resut_file = ltr_cons_path
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        log.logger.info('Start step0.3: Remove LTR redundancy')
        deredundant_for_LTR(confident_ltr_cut_path, tmp_output_dir, threads)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step0.3: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')




