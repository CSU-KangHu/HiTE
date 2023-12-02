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
        ltr_index = ltrName2ltrIndex[ltrName]
        new_name = 'ltr_' + str(ltr_index) + '-' + ltr_type
        rename_LTRs[new_name] = contigs[ltrName]
    store_fasta(rename_LTRs, rename_ltr_lib)

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
    parser.add_argument('--miu', metavar='miu',
                        help='The neutral mutation rate (per bp per ya)')


    args = parser.parse_args()
    reference = args.g
    LTR_harvest_parallel_Home = args.ltrharvest_home
    LTR_finder_parallel_Home = args.ltrfinder_home
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    recover = args.recover
    miu = args.miu

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

    resut_file = ref_rename_path + '.LTRlib.fa'
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        log.logger.info('Start step0.2: run LTR_retriever to get confident LTR')
        run_LTR_retriever(ref_rename_path, tmp_output_dir, threads, miu, log)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step0.2: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')

    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    rename_LTR(resut_file, confident_ltr_cut_path)
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




