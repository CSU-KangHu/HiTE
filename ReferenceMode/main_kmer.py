#-- coding: UTF-8 --
import argparse

import codecs
import subprocess

import datetime
import json
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count, set_start_method
#
# set_start_method('forkserver', force=True)

from module.Util import read_fasta, Logger, store_fasta, \
    get_candidate_repeats, split2cluster_normal, \
    convertToUpperCase_v1, multi_line, run_LTR_harvest, getUniqueKmer_v1, file_exist, \
    determine_repeat_boundary_v1, multi_process_TRF, multi_process_align, get_copies, flanking_copies, store_copies_v1, \
    generate_candidate_ltrs, rename_reference, store_LTR_seq_v1, store_LTR_seq, multi_process_align_and_get_copies, \
    remove_ltr_from_tir, flanking_seq, rename_fasta, run_LTR_retriever, flank_region_align_v1

#from module.judge_TIR_transposons import is_transposons

def is_transposons(filter_dup_path, reference, threads, tmp_output_dir, flanking_len, blast_program_dir, ref_index, debug, log):
    log.logger.info('determine true TIR')

    log.logger.info('------flank TIR copy and see if the flanking regions are repeated')
    starttime = time.time()
    # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
    # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
    flanking_len = 50
    similar_ratio = 0.1
    TE_type = 'tir'
    confident_copies = flank_region_align_v1(filter_dup_path, flanking_len, similar_ratio, reference, TE_type, tmp_output_dir, blast_program_dir, threads, ref_index, debug, log)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of flanking TIR copy and see if the flanking regions are repeated: %.8s s" % (dtime))

    log.logger.info('------store confident TIR sequences')
    filter_dup_names, filter_dup_contigs = read_fasta(filter_dup_path)
    if ref_index == -1:
        confident_tir_path = tmp_output_dir + '/confident_tir.rename.cons.fa'
    else:
        confident_tir_path = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
    confident_tir = {}
    for name in confident_copies.keys():
        copy_list = confident_copies[name]
        if len(copy_list) >= 2:
            confident_tir[name] = filter_dup_contigs[name]
    store_fasta(confident_tir, confident_tir_path)


if __name__ == '__main__':
    default_k_num = 31
    default_threads = int(cpu_count())
    default_fault_tolerant_bases = 200
    default_fixed_extend_base_threshold = 1000
    default_chunk_size = 200
    default_freq_threshold = 2
    default_tandem_region_cutoff = 0.5
    default_max_single_repeat_len = 30000
    default_plant = 1
    default_nested = 1
    default_recover = 0
    default_flanking_len = 50
    default_global_flanking_filter = 1
    default_debug = 0
    default_chrom_seg_length = 500000

    version_num = '1.0.1'

    # 1.parse args
    parser = argparse.ArgumentParser(description='########################## HiTE, version ' + str(version_num) + ' ##########################')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('-k', metavar='kmer size',
                        help='input kmer size, default = [ '+str(default_k_num)+' ]')
    parser.add_argument('-t', metavar='thread num',
                        help='input thread num, default = [ '+str(default_threads)+' ]')
    parser.add_argument('-a', metavar='alias name',
                        help='input alias name')
    parser.add_argument('--fault_tolerant_bases', metavar='fault_tolerant_bases',
                        help='the base number of fault tolerant in repeated kmers masking, default = [ '+str(default_fault_tolerant_bases)+' ]')
    parser.add_argument('--fixed_extend_base_threshold', metavar='fixed_extend_base_threshold',
                        help='the base number of extend base, default = [ '+str(default_fixed_extend_base_threshold)+' ]')
    parser.add_argument('--chunk_size', metavar='chunk_size',
                        help='the chunk size of large genome, default = [ '+str(default_chunk_size)+' MB ]')
    parser.add_argument('--freq_threshold', metavar='freq_threshold',
                        help='the frequency threshold of kmer, default = [ '+str(default_freq_threshold)+' ]')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='Cutoff of the raw masked repeat regarded as tandem region, default = [ '+str(default_tandem_region_cutoff)+' ]')
    parser.add_argument('--max_repeat_len', metavar='max_repeat_len',
                        help='the maximum length of repeat, default = [ ' + str(default_max_single_repeat_len) + ' ]')
    parser.add_argument('--chrom_seg_length', metavar='chrom_seg_length',
                        help='the length of genome segments, default = [ ' + str(default_chrom_seg_length) + ' ]')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='the flanking length of repeat to find the true boundary, default = [ '+str(default_flanking_len)+' ]')
    parser.add_argument('--plant', metavar='is_plant',
                        help='is it a plant genome, 1: true, 0: false. default = [ ' + str(default_plant) + ' ]')
    parser.add_argument('--remove_nested', metavar='remove_nested',
                        help='Whether to clear the nested TE, 1: true, 0: false. default = [ ' + str(default_nested) + ' ]')
    parser.add_argument('--global_flanking_filter', metavar='global_flanking_filter',
                        help='Whether to filter false positives by global flanking alignment, significantly reduce false positives '
                             'but require more memory, especially when inputting a large genome. 1: true (require more memory), 0: false. default = [ ' + str(
                            default_global_flanking_filter) + ' ]')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid repeated calculations, 1: true, 0: false. default = [ ' + str(
                            default_recover) + ' ]')
    parser.add_argument('--debug', metavar='is_debug',
                        help='Open debug mode, 1: true, 0: false. default = [ ' + str(default_debug) + ' ]')
    parser.add_argument('-o', metavar='output dir',
                        help='output dir')

    args = parser.parse_args()

    reference = args.g
    k_num = args.k
    threads = args.t
    alias = args.a
    output_dir = args.o
    fault_tolerant_bases = args.fault_tolerant_bases
    fixed_extend_base_threshold = args.fixed_extend_base_threshold
    chunk_size = args.chunk_size
    freq_threshold = args.freq_threshold
    tandem_region_cutoff = args.tandem_region_cutoff
    max_repeat_len = args.max_repeat_len
    chrom_seg_length = args.chrom_seg_length
    flanking_len = args.flanking_len
    plant = args.plant
    remove_nested = args.remove_nested
    global_flanking_filter = args.global_flanking_filter
    recover = args.recover
    debug = args.debug

    log = Logger('HiTE.log', level='debug')

    if reference is None:
        log.logger.error('\nreference path can not be empty')
        parser.print_help()
        exit(-1)
    if output_dir is None:
        output_dir = os.getcwd() + '/output'
        log.logger.warning('\noutput directory path is empty, set to: ' + str(output_dir))

    if not os.path.isabs(reference):
        reference = os.path.abspath(reference)
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    if k_num is None:
        k_num = int(default_k_num)
    else:
        k_num = int(k_num)

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    if fault_tolerant_bases is None:
        fault_tolerant_bases = default_fault_tolerant_bases
    else:
        fault_tolerant_bases = int(fault_tolerant_bases)

    if fixed_extend_base_threshold is None:
        fixed_extend_base_threshold = default_fixed_extend_base_threshold
    else:
        fixed_extend_base_threshold = int(fixed_extend_base_threshold)

    if chunk_size is None:
        chunk_size = default_chunk_size
    else:
        chunk_size = float(chunk_size)

    if freq_threshold is None:
        freq_threshold = default_freq_threshold
    else:
        freq_threshold = int(freq_threshold)

    if tandem_region_cutoff is None:
        tandem_region_cutoff = default_tandem_region_cutoff
    else:
        tandem_region_cutoff = float(tandem_region_cutoff)

    if flanking_len is None:
        flanking_len = default_flanking_len
    else:
        flanking_len = int(flanking_len)

    if plant is None:
        plant = default_plant
    else:
        plant = int(plant)

    if remove_nested is None:
        remove_nested = default_nested
    else:
        remove_nested = int(remove_nested)

    if global_flanking_filter is None:
        global_flanking_filter = default_global_flanking_filter
    else:
        global_flanking_filter = int(global_flanking_filter)

    if recover is None:
        recover = default_recover
    else:
        recover = int(recover)

    if debug is None:
        debug = default_debug
    else:
        debug = int(debug)

    if max_repeat_len is None:
        max_repeat_len = default_max_single_repeat_len
    else:
        max_repeat_len = int(max_repeat_len)

    if chrom_seg_length is None:
        chrom_seg_length = default_chrom_seg_length
    else:
        chrom_seg_length = int(chrom_seg_length)

    partitions_num = int(threads)

    # Step1. read configuration
    param_config_path = os.getcwd() + "/ParamConfig.json"
    # read param config
    with open(param_config_path, 'r') as load_f:
        param = json.load(load_f)
    load_f.close()

    i = datetime.datetime.now()
    # tmp_output_dir = output_dir + '/CRD.' + str(i.date()) + '.' + str(i.hour) + '-' + str(i.minute) + '-' + str(i.second)
    tmp_output_dir = output_dir + '/'
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    total_starttime = time.time()
    tools_dir = os.getcwd() + '/tools'

    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

    os.system('cp ' + reference + ' ' + tmp_output_dir)
    reference = tmp_output_dir + '/' + ref_filename

    # preset steps, no needs to execute time-consuming job again
    # when the job is retried due to abnormal termination.
    # 检查各个步骤的结果文件
    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    blast_program_dir = param['RMBlast_Home']
    Genome_Tools_Home = param['Genome_Tools_Home']
    LTR_retriever_Home = param['LTR_retriever_Home']
    TRF_Path = param['TRF_Path']

    # LTR_finder_parallel_Home = param['LTR_finder_parallel_Home']
    # EAHelitron = param['EAHelitron']

    LTR_finder_parallel_Home = os.getcwd() + '/bin/LTR_FINDER_parallel-master'
    EAHelitron = os.getcwd() + '/bin/EAHelitron-master'

    if blast_program_dir == '':
        (status, blast_program_path) = subprocess.getstatusoutput('which makeblastdb')
        blast_program_dir = os.path.dirname(os.path.dirname(blast_program_path))
    if Genome_Tools_Home == '':
        (status, Genome_Tools_path) = subprocess.getstatusoutput('which gt')
        Genome_Tools_Home = os.path.dirname(os.path.dirname(Genome_Tools_path))
    if LTR_retriever_Home == '':
        (status, LTR_retriever_path) = subprocess.getstatusoutput('which LTR_retriever')
        LTR_retriever_Home = os.path.dirname(LTR_retriever_path)
    if TRF_Path == '':
        (status, TRF_Path) = subprocess.getstatusoutput('which trf')


    log.logger.info('\n-------------------------------------------------------------------------------------------\n'
                    'HiTE, version ' + str(version_num) + '\n'
                    'Copyright (C) 2022 Kang Hu ( kanghu@csu.edu.cn )\n'
                    'Hunan Provincial Key Lab on Bioinformatics, School of Computer Science and \n'
                    'Engineering, Central South University, Changsha 410083, P.R. China.\n'
                    '-------------------------------------------------------------------------------------------')

    log.logger.info('\nParameters configuration\n'
                    '====================================System settings========================================\n'
                    '  [Setting] Reference sequences / assemblies path = [ ' + str(reference) + ' ]\n'
                    '  [Setting] Alias = [ ' + str(alias) + ' ]\n'
                    '  [Setting] The K-mer Size = [ ' + str(k_num) + 'bp]  Default( ' + str(default_k_num) + ' )\n'
                    '  [Setting] Threads = [ ' + str(threads) + ' ]  Default( ' + str(default_threads) + ' )\n'
                    '  [Setting] Fault tolerant bases threshold = [ ' + str(fault_tolerant_bases) + ' ] Default( ' + str(default_fault_tolerant_bases) + ' )\n'
                    '  [Setting] Fixed extend bases threshold = [ ' + str(fixed_extend_base_threshold) + ' ] Default( ' + str(default_fixed_extend_base_threshold) + ' )\n'
                    
                    '  [Setting] The chunk size of large genome = [ ' + str(chunk_size) + ' ] MB Default( ' + str(default_chunk_size) + ' ) MB\n'
                    '  [Setting] The frequency threshold of kmer = [ ' + str(freq_threshold) + ' ] Default( ' + str(default_freq_threshold) + ' )\n'
                    '  [Setting] Cutoff of the repeat regarded as tandem sequence = [ ' + str(tandem_region_cutoff) + ' ] Default( ' + str(default_tandem_region_cutoff) + ' )\n'
                    '  [Setting] Maximum length of TE = [ ' + str(max_repeat_len) + ' ]  Default( ' + str(default_max_single_repeat_len) + ' )\n'
                    '  [Setting] The length of genome segments = [ ' + str(chrom_seg_length) + ' ]  Default( ' + str(default_chrom_seg_length) + ' )\n'
                    '  [Setting] Flanking length of TE = [ ' + str(flanking_len) + ' ]  Default( ' + str(default_flanking_len) + ' )\n'
                    '  [Setting] Is plant genome = [ ' + str(plant) + ' ]  Default( ' + str(default_plant) + ' )\n'
                    '  [Setting] Remove nested = [ ' + str(remove_nested) + ' ]  Default( ' + str(default_nested) + ' )\n'
                    '  [Setting] Global flanking filter = [ ' + str(global_flanking_filter) + ' ]  Default( ' + str(default_global_flanking_filter) + ' )\n'
                    '  [Setting] recover = [ ' + str(recover) + ' ]  Default( ' + str(default_recover) + ' )\n'
                    '  [Setting] debug = [ ' + str(debug) + ' ]  Default( ' + str(default_debug) + ' )\n'
                    '  [Setting] Output Directory = [' + str(output_dir) + ']'

                    '  [Setting] Blast Program Home = [' + str(blast_program_dir) + ']\n'
                    '  [Setting] Genome Tools Program Home = [' + str(Genome_Tools_Home) + ']\n'
                    '  [Setting] LTR_retriever Program Home = [' + str(LTR_retriever_Home) + ']\n'
                    '  [Setting] TRF Program Path = [' + str(TRF_Path) + ']'
                    )


    TRsearch_dir = tools_dir
    test_home = os.getcwd() + '/module'

    pipeline_starttime = time.time()
    # 我们将大的基因组划分成多个小的基因组，每个小基因组500M，分割来处理
    starttime = time.time()

    # --------------------------------------------------------------------------------------
    # Step1. dsk get unique kmers, whose frequency >= 2
    log.logger.info('Start Splitting Reference into chunks')
    # using multiple threads to gain speed
    reference_pre = convertToUpperCase_v1(reference)
    reference_tmp = multi_line(reference_pre, chrom_seg_length, k_num)
    cut_references = []
    cur_ref_contigs = {}
    cur_base_num = 0
    ref_index = 0
    with open(reference_tmp, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ref_name = parts[0].replace('>', '')
            start = parts[1]
            seq = parts[2]
            new_ref_name = ref_name + '$' + start
            cur_ref_contigs[new_ref_name] = seq
            cur_base_num += len(line)
            if cur_base_num >= chunk_size * 1024 * 1024:
                # store references
                cur_ref_path = reference + '.cut'+str(ref_index)+'.fa'
                store_fasta(cur_ref_contigs, cur_ref_path)
                cut_references.append(cur_ref_path)
                cur_ref_contigs = {}
                cur_base_num = 0
                ref_index += 1
        if len(cur_ref_contigs) > 0:
            cur_ref_path = reference + '.cut' + str(ref_index) + '.fa'
            store_fasta(cur_ref_contigs, cur_ref_path)
            cut_references.append(cur_ref_path)

    for ref_index, cut_reference in enumerate(cut_references):
        (cut_ref_dir, cut_ref_filename) = os.path.split(cut_reference)
        (cut_ref_name, cut_ref_extension) = os.path.splitext(cut_ref_filename)
        log.logger.info('Round of chunk: ' + str(ref_index))

        log.logger.info('Start step1: K-mer Based De Novo TE Searching')
        connected_repeats = {}
        resut_file = tmp_output_dir + '/repeats_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            log.logger.info('Start step1.1: Running RepeatMasking Algorithm')
            starttime = time.time()
            dsk_h5_path = cut_ref_name + '.h5'
            unique_kmer_path = tmp_output_dir + '/kmer.txt'
            dsk_cmd1 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk -file ' + cut_reference + ' -kmer-size ' + str(k_num) + ' -abundance-min ' + str(freq_threshold)
            dsk_cmd2 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk2ascii -file ' + dsk_h5_path + ' -out ' + unique_kmer_path
            #log.logger.debug(dsk_cmd1)
            os.system(dsk_cmd1 + ' > /dev/null 2>&1')
            #log.logger.debug(dsk_cmd2)
            os.system(dsk_cmd2 + ' > /dev/null 2>&1')

            # --------------------------------------------------------------------------------------
            # Step2. each thread process a batch of kmers

            # create thread pool, use multiple processes to execute
            kmer_segments = []
            with open(unique_kmer_path, 'r') as f_r:
                for line in f_r:
                    line = line.replace('\n', '')
                    kmer_segments.append(line)
            kmer_segments_cluster = split2cluster_normal(kmer_segments, partitions_num)

            ex = ProcessPoolExecutor(partitions_num)
            objs = []
            for partiton_index in kmer_segments_cluster.keys():
                cur_kmer_segments = kmer_segments_cluster[partiton_index]
                obj = ex.submit(getUniqueKmer_v1, cur_kmer_segments, partiton_index)
                objs.append(obj)
            ex.shutdown(wait=True)
            unique_kmer_map = {}
            for obj in as_completed(objs):
                for kmer in obj.result():
                    unique_kmer_map[kmer] = 1

            # --------------------------------------------------------------------------------------
            # Step3. get candidate repeats
            cur_connected_repeats = get_candidate_repeats(cut_reference, k_num, partitions_num, unique_kmer_map, fault_tolerant_bases, tmp_output_dir, log)

            for name in cur_connected_repeats.keys():
                repeat_list = cur_connected_repeats[name]
                connected_repeats[name+'-'+str(ref_index)] = repeat_list

            # generate repeats.fa and connected_regions
            repeats_path = tmp_output_dir + '/repeats_'+str(ref_index)+'.fa'
            node_index = 0
            with open(repeats_path, 'w') as f_save:
                for ref_name in connected_repeats.keys():
                    repeat_list = connected_repeats[ref_name]
                    for repeat_item in repeat_list:
                        start_pos = repeat_item[0]
                        end_pos = repeat_item[1]
                        # 如果序列过长，就将序列进行切割。过长通常是由于k_num取的过小导致
                        # 按照100k进行切割
                        repeat = repeat_item[2]
                        long_len = 100*1000
                        if len(repeat) >= long_len:
                            new_start_pos = start_pos
                            while len(repeat) >= long_len:
                                new_end_pos = new_start_pos + long_len
                                new_query_name = 'P_' + str(ref_index) + '-N_' + str(node_index) + '-s_' + str(ref_name) + '-' + str(new_start_pos) + '-' + str(new_end_pos)
                                new_repeat = repeat[0: long_len]
                                if len(new_repeat) >= 80:
                                    f_save.write('>' + new_query_name + '\n' + new_repeat + '\n')
                                new_start_pos = new_end_pos + 1
                                repeat = repeat[long_len:]
                            if len(repeat) >= 80:
                                new_query_name = 'P_' + str(ref_index) + '-N_' + str(node_index) + '-s_' + str(ref_name) + '-' + str(new_start_pos) + '-' + str(end_pos)
                                f_save.write('>' + new_query_name + '\n' + repeat + '\n')
                        else:
                            query_name = 'P_' + str(ref_index) + '-N_' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
                            repeat = repeat_item[2]
                            if len(repeat) >= 80:
                                f_save.write('>' + query_name + '\n' + repeat + '\n')
                                node_index += 1

            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of RepeatMasking Algorithm: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')


        #longest_repeats_copies_file = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.copies'
        longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
        longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
        resut_file = tmp_output_dir + '/longest_repeats_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file) or not file_exist(longest_repeats_flanked_path):
            starttime = time.time()
            log.logger.info('Start 2.2: Coarse-grained boundary mapping')
            log.logger.info('------generate longest_repeats.fa')
            repeats_path = tmp_output_dir + '/repeats_'+str(ref_index)+'.fa'
            # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
            determine_repeat_boundary_v1(repeats_path, longest_repeats_path, blast_program_dir,
                                         fixed_extend_base_threshold, max_repeat_len, reference, tmp_output_dir, debug, threads)

            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of generating longest_repeats.fa: %.8s s" % (dtime))

            starttime = time.time()
            trf_dir = tmp_output_dir + '/trf_temp'
            (repeat_dir, repeat_filename) = os.path.split(longest_repeats_path)
            (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
            repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.filter_tandem.fa'
            multi_process_TRF(longest_repeats_path, repeats_path, TRF_Path, trf_dir, tandem_region_cutoff,
                              threads=threads)
            if debug == 0:
                # remove temp dir
                os.system('rm -rf ' + trf_dir)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of filtering tandem repeat in longest_repeats.fa: %.8s s" % (dtime))

            longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.filter_tandem.fa'
            longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
            flanking_seq(longest_repeats_path, longest_repeats_flanked_path, reference, flanking_len)


            # longest_repeats_path_rename = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.rename.fa'
            # contignames, contigs = read_fasta(longest_repeats_path)
            # node_index = 0
            # with open(longest_repeats_path_rename, 'w') as f_save:
            #     for name in contignames:
            #         seq = contigs[name]
            #         f_save.write('>N_' + str(node_index) + '\n' + seq + '\n')
            #         node_index += 1
            #
            # log.logger.info('------get copies of longest_repeats_' + str(ref_index) + '.rename.fa')
            # starttime = time.time()
            # # 1.因为现在的TIR序列不能准确的确定边界，因此下面的步骤都是为了能够准确获得边界
            # valid_TIR_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.rename.fa'
            # blastnResults_path = tmp_output_dir + '/valid_tir.ref.out'
            # valid_tir_blast_dir = tmp_output_dir + '/valid_tir_blast'
            # all_copies = multi_process_align_and_get_copies(valid_TIR_path, cut_reference, blast_program_dir, valid_tir_blast_dir, threads)
            # # multi_process_align(valid_TIR_path, cut_reference, blastnResults_path, blast_program_dir, valid_tir_blast_dir,
            # #                     threads)
            # # all_copies = get_copies(blastnResults_path, valid_TIR_path, cut_reference, query_coverage=0.99, threads=threads)
            # endtime = time.time()
            # dtime = endtime - starttime
            # log.logger.info("Running time of getting copies: %.8s s" % (dtime))
            #
            # log.logger.info('------flanking copies')
            # starttime = time.time()
            # # 需要多个拷贝来判断是否是真实的TSD
            # flanked_all_copies = flanking_copies(all_copies, valid_TIR_path, cut_reference, flanking_len, copy_num=20, query_coverage=0.99)
            # endtime = time.time()
            # dtime = endtime - starttime
            # log.logger.info("Running time of flanking copies: %.8s s" % (dtime))

            # # 存储all copies，在后面的judge TIR 和 Helitron会用到
            # with codecs.open(longest_repeats_flanked_copies_file, 'w', encoding='utf-8') as f:
            #     json.dump(flanked_all_copies, f)
            #
            # copy_info_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.rename.fa.copies.info'
            # store_copies_v1(flanked_all_copies, copy_info_path)

        else:
            log.logger.info(resut_file + ' exists, skip...')

        # raw_candidate_ltrs_file = tmp_output_dir + '/raw_candidate_ltr_' + str(ref_index) + '.fa'
        # resut_file = raw_candidate_ltrs_file
        # if not is_recover or not file_exist(resut_file):
        #     log.logger.info('------generate candidate ltr sequences')
        #     starttime = time.time()
        #
        #     file = open(longest_repeats_copies_file, 'r')
        #     js = file.read()
        #     all_copies = json.loads(js)
        #
        #     raw_candidate_ltrs = generate_candidate_ltrs(all_copies, reference, flanking_len)
        #     endtime = time.time()
        #     dtime = endtime - starttime
        #     log.logger.info("Running time of generating candidate ltr sequences: %.8s s" % (dtime))
        #     # raw_candidate_ltrs，在后面的judge LTR 会用到
        #     with codecs.open(raw_candidate_ltrs_file, 'w', encoding='utf-8') as f:
        #         json.dump(raw_candidate_ltrs, f)


        # resut_file = tmp_output_dir + '/new_confident_ltr_' + str(ref_index) + '.fa'
        # if not is_recover or not file_exist(resut_file):
        #     starttime = time.time()
        #     log.logger.info('Start step5: determine fine-grained new LTR')
        #     # 识别LTR转座子
        #     ltr_identification_command = 'cd ' + test_home + ' && python ' + test_home + '/search_LTR_transposons.py -g ' \
        #                                  + cut_reference + ' --raw_candidate_ltrs_file ' + raw_candidate_ltrs_file \
        #                                  + ' -t ' + str(threads) + ' --TRsearch_dir ' + TRsearch_dir \
        #                                  + ' --tmp_output_dir ' + tmp_output_dir + ' --blast_program_dir ' \
        #                                  + blast_program_dir + ' --ref_index ' + str(ref_index) \
        #                                  + ' --plant ' + str(plant) + ' --flanking_len ' + str(flanking_len)
        #     os.system(ltr_identification_command)
        #     endtime = time.time()
        #     dtime = endtime - starttime
        #     log.logger.info("Running time of step5: %.8s s" % (dtime))
        # else:
        #     log.logger.info(resut_file + ' exists, skip...')

        longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
        resut_file = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1.3: determine fine-grained TIR')
            # 识别TIR转座子
            tir_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_TIR_transposons.py -g ' \
                                         + cut_reference + ' --seqs ' + longest_repeats_flanked_path\
                                         + ' -t ' + str(threads)+' --TRsearch_dir ' + TRsearch_dir \
                                         + ' --tmp_output_dir ' + tmp_output_dir + ' --blast_program_dir ' \
                                         + blast_program_dir + ' --TRF_Path ' + TRF_Path + ' --tandem_region_cutoff ' \
                                         + str(tandem_region_cutoff) + ' --ref_index ' + str(ref_index) \
                                         + ' --plant ' + str(plant) + ' --flanking_len ' + str(flanking_len) \
                                         + ' --debug ' + str(debug)
            os.system(tir_identification_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1.3: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        resut_file = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1.4: determine fine-grained Helitron')
            # 识别Helitron转座子
            helitron_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Helitron_transposons.py --seqs ' \
                                              + longest_repeats_flanked_path + ' -g ' + cut_reference + ' -t ' + str(threads) \
                                              + ' --tmp_output_dir ' + tmp_output_dir + ' --EAHelitron ' + EAHelitron \
                                              + ' --blast_program_dir ' + blast_program_dir \
                                              + ' --ref_index ' + str(ref_index) + ' --flanking_len ' + str(flanking_len) \
                                              + ' --debug ' + str(debug)

            # HSDIR = '/public/home/hpc194701009/repeat_detect_tools/TrainingSet'
            # HSJAR = '/public/home/hpc194701009/repeat_detect_tools/HelitronScanner/HelitronScanner.jar'
            # helitron_identification_command = 'cd ' + test_home + ' && python ' + test_home + '/judge_Helitron_transposons.py --copies ' \
            #                              + longest_repeats_copies_file + ' -g ' + cut_reference + ' -t ' + str(threads)+' --HSDIR ' + HSDIR \
            #                              + ' --HSJAR ' + HSJAR + ' --tmp_output_dir ' + tmp_output_dir \
            #                              + ' --blast_program_dir ' + blast_program_dir \
            #                              + ' --ref_index ' + str(ref_index)
            os.system(helitron_identification_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1.4: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        resut_file = tmp_output_dir + '/confident_other_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1.5: homology-based other TE searching')
            #同源搜索其他转座子
            other_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Other_transposons.py ' \
                                           + ' --seqs ' + longest_repeats_flanked_path\
                                           + ' -t ' + str(threads) + ' --blast_program_dir ' + blast_program_dir \
                                           + ' --tmp_output_dir ' + tmp_output_dir + ' --query_coverage ' + str(0.8) \
                                           + ' --subject_coverage ' + str(0) + ' --ref_index ' + str(ref_index)
            os.system(other_identification_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1.5: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        #合并Helitron、TIR、other转座子
        cur_confident_TE_path = tmp_output_dir + '/confident_TE_'+str(ref_index)+'.fa'
        cur_confident_helitron_path = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
        cur_confident_other_path = tmp_output_dir + '/confident_other_'+str(ref_index)+'.fa'
        os.system('cat ' + cur_confident_helitron_path + ' > ' + cur_confident_TE_path)
        os.system('cat ' + cur_confident_other_path + ' >> ' + cur_confident_TE_path)

    # 过滤TIR候选序列中的LTR转座子（intact LTR or LTR terminals or LTR internals）
    # 1.1 合并所有parts的TIR序列
    confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    confident_other_path = tmp_output_dir + '/confident_other.fa'
    os.system('rm -f ' + confident_tir_path)
    os.system('rm -f ' + confident_other_path)
    for ref_index, ref_rename_path in enumerate(cut_references):
        cur_confident_tir_path = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
        cur_confident_other_path = tmp_output_dir + '/confident_other_' + str(ref_index) + '.fa'
        os.system('cat ' + cur_confident_tir_path + ' >> ' + confident_tir_path)
        os.system('cat ' + cur_confident_other_path + ' >> ' + confident_other_path)

    log.logger.info('Start step2: Structural Based LTR Searching')
    # 1.重命名reference文件
    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

    ref_rename_path = tmp_output_dir + '/' + ref_name + '.rename.fa'
    rename_reference(reference, ref_rename_path)

    backjob = None
    resut_file = tmp_output_dir + '/genome_all.fa.harvest.scn'
    if not is_recover or not file_exist(resut_file):
        log.logger.info('Start step2.1: Running LTR_harvest')
        # -------------------------------Stage01: this stage is used to generate kmer coverage repeats-------------------------------
        # --------------------------------------------------------------------------------------
        # run LTRharvest background job
        backjob = multiprocessing.Process(target=run_LTR_harvest,
                                          args=(Genome_Tools_Home, ref_rename_path, tmp_output_dir, log,))
        backjob.start()
    else:
        log.logger.info(resut_file + ' exists, skip...')

    resut_file = ref_rename_path + '.finder.combine.scn'
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        log.logger.info('Start step2.2: Running LTR finder parallel to obtain candidate LTRs')

        # 运行LTR_finder_parallel来获取候选的LTR序列
        # 2.运行LTR_finder_parallel
        LTR_finder_parallel_command = 'perl ' + LTR_finder_parallel_Home +'/LTR_FINDER_parallel -harvest_out -seq ' + ref_rename_path + ' -threads ' + str(
            threads)
        os.system('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')

        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of obtaining candidate LTRs: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')

    # 合并LTR_harvest+LTR_finder结果，输入到LTR_retriever
    if backjob is not None:
        backjob.join()
    ltrharvest_output = tmp_output_dir + '/genome_all.fa.harvest.scn'
    ltrfinder_output = ref_rename_path + '.finder.combine.scn'
    ltr_output = tmp_output_dir + '/genome_all.fa.rawLTR.scn'
    os.system('cat ' + ltrharvest_output + ' ' + ltrfinder_output + ' > ' + ltr_output)

    resut_file = ref_rename_path + '.LTRlib.fa'
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        log.logger.info('Start step2.3: run LTR_retriever to get confident LTR')
        run_LTR_retriever(LTR_retriever_Home, ref_rename_path, tmp_output_dir, threads, log)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step2.3: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')

    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    os.system('cp ' + resut_file + ' ' + confident_ltr_cut_path)


    # 1.2 confident_ltr_cut_path比对到TIR候选序列上，并且过滤掉出现在LTR库中的TIR序列
    temp_dir = tmp_output_dir + '/tir_blast_ltr'
    all_copies = multi_process_align_and_get_copies(confident_ltr_cut_path, confident_tir_path, blast_program_dir, temp_dir, 'tir', threads, query_coverage=0.8)
    remove_ltr_from_tir(confident_ltr_cut_path, confident_tir_path, all_copies)

    # # 1.3 比对到confident_ltr上，并且过滤掉出现在LTR库中的TIR序列
    # temp_dir = tmp_output_dir + '/tir_blast_ltr'
    # all_copies = multi_process_align_and_get_copies(confident_ltr_path, confident_tir_path, blast_program_dir,
    #                                                 temp_dir, 'tir', threads, query_coverage=0.8)
    # remove_ltr_from_tir(confident_ltr_path, confident_tir_path, all_copies)

    # 1.4 生成一致性tir序列
    confident_tir_rename_path = tmp_output_dir + '/confident_tir.rename.fa'
    rename_fasta(confident_tir_path, confident_tir_rename_path)

    confident_tir_rename_consensus = tmp_output_dir + '/confident_tir.rename.cons.fa'
    tools_dir = os.getcwd() + '/tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_tir_rename_path + ' -o ' + confident_tir_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    # 如果切分成了多个块，TIR需要重新flank_region_align_v1到整个基因组，以过滤掉那些在分块中未能过滤掉的假阳性。
    if len(cut_references) > 1 and global_flanking_filter == 1:
        ref_index = -1
        is_transposons(confident_tir_rename_consensus, reference, threads, tmp_output_dir, flanking_len, blast_program_dir,
                       ref_index, debug, log)

    # 1.5 解开TIR中包含的nested TE
    clean_tir_path = tmp_output_dir + '/confident_tir.clean.fa'
    remove_nested_command = 'cd ' + test_home + ' && python3 remove_nested_lib.py ' \
                            + ' -t ' + str(threads) + ' --blast_program_dir ' + blast_program_dir \
                            + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(5) \
                            + ' --input1 ' + confident_tir_rename_consensus \
                            + ' --input2 ' + confident_tir_rename_consensus \
                            + ' --output ' + clean_tir_path
    os.system(remove_nested_command)

    # # cd-hit -aS 0.95 -c 0.8合并一些冗余序列
    # clean_tir_consensus = tmp_output_dir + '/confident_tir.clean.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + clean_tir_path + ' -o ' + clean_tir_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)

    # 1.6 生成一致性other序列
    confident_other_rename_path = tmp_output_dir + '/confident_other.rename.fa'
    rename_fasta(confident_other_path, confident_other_rename_path)

    confident_other_rename_consensus = tmp_output_dir + '/confident_other.rename.cons.fa'
    tools_dir = os.getcwd() + '/tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_other_rename_path + ' -o ' + confident_other_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    # 合并所有parts的TE、TIR转座子
    confident_TE_path = tmp_output_dir + '/confident_TE.fa'
    os.system('rm -f ' + confident_TE_path)
    for ref_index, cut_reference in enumerate(cut_references):
        cur_confident_TE_path = tmp_output_dir + '/confident_TE_' + str(ref_index) + '.fa'
        os.system('cat ' + cur_confident_TE_path + ' >> ' + confident_TE_path)
    os.system('cat ' + clean_tir_path + ' >> ' + confident_TE_path)

    # 解开LTR内部包含的nested TE，然后把解开后的LTR合并到TE库中
    # 获取LTR的内部序列
    confident_ltr_terminal_path = tmp_output_dir + '/confident_ltr_cut.terminal.fa'
    confident_ltr_internal_path = tmp_output_dir + '/confident_ltr_cut.internal.fa'
    ltr_names, ltr_contigs = read_fasta(confident_ltr_cut_path)
    ltr_internal_contigs = {}
    ltr_terminal_contigs = {}
    for name in ltr_names:
        if name.__contains__('_INT#'):
            ltr_internal_contigs[name] = ltr_contigs[name]
        else:
            ltr_terminal_contigs[name] = ltr_contigs[name]
    store_fasta(ltr_internal_contigs, confident_ltr_internal_path)
    store_fasta(ltr_terminal_contigs, confident_ltr_terminal_path)

    clean_ltr_internal_path = tmp_output_dir + '/confident_ltr_cut.internal.clean.fa'
    if remove_nested == 1:
        starttime = time.time()
        log.logger.info('Start step2.4: remove nested TE in LTR internal')
        # 将所有的LTR序列暂时加入到TE序列中，用来解开nested TE
        temp_confident_TE_path = tmp_output_dir + '/confident_TE.temp.fa'
        os.system('cat ' + confident_ltr_cut_path + ' > ' + temp_confident_TE_path)
        os.system('cat ' + confident_TE_path + ' >> ' + temp_confident_TE_path)
        remove_nested_command = 'cd ' + test_home + ' && python3 remove_nested_lib.py ' \
                                + ' -t ' + str(threads) + ' --blast_program_dir ' + blast_program_dir \
                                + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(5) \
                                + ' --input1 ' + temp_confident_TE_path \
                                + ' --input2 ' + confident_ltr_internal_path \
                                + ' --output ' + clean_ltr_internal_path
        os.system(remove_nested_command)

        os.system('cat ' + confident_ltr_terminal_path + ' > ' + confident_ltr_cut_path)
        os.system('cat ' + clean_ltr_internal_path + ' >> ' + confident_ltr_cut_path)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step2.4: %.8s s" % (dtime))
    os.system('cat ' + confident_ltr_cut_path + ' >> ' + confident_TE_path)

    confident_ltr_cut_consensus = tmp_output_dir + '/confident_ltr_cut.cons.fa'
    tools_dir = os.getcwd() + '/tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_ltr_cut_path + ' -o ' + confident_ltr_cut_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    starttime = time.time()
    log.logger.info('Start step3: generate non-redundant library')
    generate_lib_command = 'cd ' + test_home + ' && python3 ' + test_home + '/get_nonRedundant_lib.py' \
                           + ' -t ' + str(threads) + ' --tmp_output_dir ' + tmp_output_dir \
                           + ' --sample_name ' + alias + ' --blast_program_dir ' + blast_program_dir
    os.system(generate_lib_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of step3: %.8s s" % (dtime))


    pipeline_endtime = time.time()
    dtime = pipeline_endtime - pipeline_starttime
    log.logger.info("Running time of the whole pipeline: %.8s s" % (dtime))
