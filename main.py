#-- coding: UTF-8 --
import argparse

import codecs
import re
import subprocess

import datetime
import json
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

from module.Util import read_fasta, Logger, store_fasta, \
    get_candidate_repeats, split2cluster_normal, \
    convertToUpperCase_v1, multi_line, run_LTR_harvest, getUniqueKmer_v1, file_exist, \
    determine_repeat_boundary_v1, multi_process_TRF, multi_process_align, get_copies, flanking_copies, store_copies_v1, \
    generate_candidate_ltrs, rename_reference, store_LTR_seq_v1, store_LTR_seq, multi_process_align_and_get_copies, \
    remove_ltr_from_tir, flanking_seq, rename_fasta, run_LTR_retriever, flank_region_align_v1, \
    determine_repeat_boundary_v2, determine_repeat_boundary_v3


#from module.judge_TIR_transposons import is_transposons

if __name__ == '__main__':
    default_threads = int(cpu_count())
    default_fixed_extend_base_threshold = 1000
    default_chunk_size = 400
    default_tandem_region_cutoff = 0.5
    default_max_single_repeat_len = 30000
    default_plant = 1
    default_recover = 0
    default_flanking_len = 50
    default_debug = 0
    default_chrom_seg_length = 500000
    default_classified = 1
    default_miu = str(1.3e-8)

    version_num = '2.0.3'
 
    # 1.parse args
    parser = argparse.ArgumentParser(description='########################## HiTE, version ' + str(version_num) + ' ##########################')
    parser.add_argument('--genome', metavar='genome', help='Input genome assembly path')
    parser.add_argument('--thread', metavar='thread_num', help='Input thread num, default = [ '+str(default_threads)+' ]')
    parser.add_argument('--chunk_size', metavar='chunk_size', help='The chunk size of large genome, default = [ ' + str(default_chunk_size) + ' MB ]')
    parser.add_argument('--miu', metavar='miu', help='The neutral mutation rate (per bp per ya), default = [ ' + str(default_miu) + ' MB ]')
    parser.add_argument('--plant', metavar='is_plant', help='Is it a plant genome, 1: true, 0: false. default = [ ' + str(default_plant) + ' ]')
    parser.add_argument('--classified', metavar='is_classified', help='Whether to classify TE models, HiTE uses RepeatClassifier from RepeatModeler to classify TEs, 1: true, 0: false. default = [ ' + str(default_classified) + ' ]')
    parser.add_argument('--recover', metavar='is_recover', help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ ' + str(default_recover) + ' ]')
    parser.add_argument('--debug', metavar='is_debug', help='Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ ' + str(default_debug) + ' ]')
    parser.add_argument('--outdir', metavar='output_dir', help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')

    parser.add_argument('--flanking_len', metavar='flanking_len', help='The flanking length of candidates to find the true boundaries, default = [ ' + str(default_flanking_len) + ' ]')
    parser.add_argument('--fixed_extend_base_threshold', metavar='fixed_extend_base_threshold', help='The length of variation can be tolerated during pairwise alignment, default = [ '+str(default_fixed_extend_base_threshold)+' ]')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff', help='Cutoff of the candidates regarded as tandem region, default = [ '+str(default_tandem_region_cutoff)+' ]')
    parser.add_argument('--max_repeat_len', metavar='max_repeat_len', help='The maximum length of a single repeat, default = [ ' + str(default_max_single_repeat_len) + ' ]')
    parser.add_argument('--chrom_seg_length', metavar='chrom_seg_length', help='The length of genome segments, default = [ ' + str(default_chrom_seg_length) + ' ]')

    args = parser.parse_args()

    reference = args.genome
    threads = args.thread
    # alias = args.a
    output_dir = args.outdir
    fixed_extend_base_threshold = args.fixed_extend_base_threshold
    chunk_size = args.chunk_size
    tandem_region_cutoff = args.tandem_region_cutoff
    max_repeat_len = args.max_repeat_len
    chrom_seg_length = args.chrom_seg_length
    flanking_len = args.flanking_len
    plant = args.plant
    classified = args.classified
    miu = args.miu
    recover = args.recover
    debug = args.debug

    i = datetime.datetime.now()
    # tmp_output_dir = output_dir + '/CRD.' + str(i.date()) + '.' + str(i.hour) + '-' + str(i.minute) + '-' + str(i.second)
    tmp_output_dir = os.path.abspath(output_dir + '/')
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')

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

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    if fixed_extend_base_threshold is None:
        fixed_extend_base_threshold = default_fixed_extend_base_threshold
    else:
        fixed_extend_base_threshold = int(fixed_extend_base_threshold)

    if chunk_size is None:
        chunk_size = default_chunk_size
    else:
        chunk_size = float(chunk_size)

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

    if classified is None:
        classified = default_classified
    else:
        classified = int(classified)
    
    if miu is None:
        miu = default_miu
    else:
        miu = str(miu)

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

    total_starttime = time.time()
    tools_dir = os.getcwd() + '/tools'

    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

    os.system('cp ' + reference + ' ' + tmp_output_dir)
    reference = tmp_output_dir + '/' + ref_filename

    # preset steps, no needs to execute time-consuming job again
    # when the job is retried due to abnormal termination.
    # ?????????????????????????????????
    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    blast_program_dir = param['RMBlast_Home']
    Genome_Tools_Home = param['Genome_Tools_Home']
    LTR_retriever_Home = param['LTR_retriever_Home']
    RepeatModeler_Home = param['RepeatModeler_Home']


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
    if RepeatModeler_Home == '':
        (status, RepeatClassifier_path) = subprocess.getstatusoutput('which RepeatClassifier')
        RepeatModeler_Home = os.path.dirname(RepeatClassifier_path)

    if blast_program_dir == '' or Genome_Tools_Home == '' or LTR_retriever_Home == '':
        print('Error configuration: please check the "ParamConfig.json" file. '
              'You should not meet this error if you install HiTE with Conda; '
              'if you do, please check whether you have installed all the packages required in "environment.yml" file.')
        sys.exit(-1)

    if RepeatModeler_Home == '' and classified == 1:
        print('Error configuration: You have not configured RepeatModeler2, so please run HiTE with "--classified 0"')
        sys.exit(-1)


    log.logger.info('\n-------------------------------------------------------------------------------------------\n'
                    'HiTE, version ' + str(version_num) + '\n'
                    'Copyright (C) 2022 Kang Hu ( kanghu@csu.edu.cn )\n'
                    'Hunan Provincial Key Lab on Bioinformatics, School of Computer Science and \n'
                    'Engineering, Central South University, Changsha 410083, P.R. China.\n'
                    '-------------------------------------------------------------------------------------------')

    log.logger.info('\nParameters configuration\n'
                    '====================================System settings========================================\n'
                    '  [Setting] Reference sequences / assemblies path = [ ' + str(reference) + ' ]\n'
                    '  [Setting] Is classified = [ ' + str(classified) + ' ] Default( ' + str(default_classified) + ' )\n'
                    '  [Setting] The neutral mutation rate (per bp per ya) = [ ' + str(miu) + ' ] Default( ' + str(default_miu) + ' )\n'
                    '  [Setting] Threads = [ ' + str(threads) + ' ]  Default( ' + str(default_threads) + ' )\n'
                    '  [Setting] The chunk size of large genome = [ ' + str(chunk_size) + ' ] MB Default( ' + str(default_chunk_size) + ' ) MB\n'
                    '  [Setting] Is plant genome = [ ' + str(plant) + ' ]  Default( ' + str(default_plant) + ' )\n'
                    '  [Setting] recover = [ ' + str(recover) + ' ]  Default( ' + str(default_recover) + ' )\n'
                    '  [Setting] debug = [ ' + str(debug) + ' ]  Default( ' + str(default_debug) + ' )\n'
                    '  [Setting] Output Directory = [' + str(output_dir) + ']'
                                                                                                                                                                                                           
                    '  [Setting] Fixed extend bases threshold = [ ' + str(fixed_extend_base_threshold) + ' ] Default( ' + str(default_fixed_extend_base_threshold) + ' )\n'
                    '  [Setting] Flanking length of TE = [ ' + str(flanking_len) + ' ]  Default( ' + str(default_flanking_len) + ' )\n'
                    '  [Setting] Cutoff of the repeat regarded as tandem sequence = [ ' + str(tandem_region_cutoff) + ' ] Default( ' + str(default_tandem_region_cutoff) + ' )\n'
                    '  [Setting] Maximum length of TE = [ ' + str(max_repeat_len) + ' ]  Default( ' + str(default_max_single_repeat_len) + ' )\n'
                    '  [Setting] The length of genome segments = [ ' + str(chrom_seg_length) + ' ]  Default( ' + str(default_chrom_seg_length) + ' )\n'
                                                                                                                                                 
                    '  [Setting] Blast Program Home = [' + str(blast_program_dir) + ']\n'
                    '  [Setting] Genome Tools Program Home = [' + str(Genome_Tools_Home) + ']\n'
                    '  [Setting] LTR_retriever Program Home = [' + str(LTR_retriever_Home) + ']\n'
                    '  [Setting] RepeatModeler Program Home = [' + str(RepeatModeler_Home) + ']'
                    )


    TRsearch_dir = tools_dir
    test_home = os.getcwd() + '/module'
    library_dir = os.getcwd() + '/library'

    pipeline_starttime = time.time()
    # ???????????????????????????????????????????????????????????????????????????500M??????????????????

    # --------------------------------------------------------------------------------------
    starttime = time.time()
    log.logger.info('Start step0: Splitting genome assembly into chunks')
    # ??????TIR?????????
    split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
                                 + reference + ' --tmp_output_dir ' + tmp_output_dir \
                                 + ' --chrom_seg_length ' + str(chrom_seg_length) + ' --chunk_size ' + str(chunk_size)
    os.system(split_genome_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of step0: %.8s s" % (dtime))

    reg_str = 'genome.cut(\d).fa$'
    cut_references = []
    for filename in os.listdir(tmp_output_dir):
        if re.match(reg_str, filename) is not None:
            cut_references.append(tmp_output_dir + '/' + filename)

    for ref_index, cut_reference in enumerate(cut_references):
        log.logger.info('Round of chunk: ' + str(ref_index))

        longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
        longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
        resut_file = longest_repeats_path
        if not is_recover or not file_exist(resut_file) or not file_exist(longest_repeats_flanked_path):
            starttime = time.time()
            log.logger.info('Start 1.1: Coarse-grained boundary mapping')
            coarse_boundary_command = 'cd ' + test_home + ' && python3 ' + test_home + '/coarse_boundary.py ' \
                                   + '-g ' + cut_reference + ' --tmp_output_dir ' + tmp_output_dir \
                                   + ' --fixed_extend_base_threshold ' + str(fixed_extend_base_threshold) \
                                   + ' --max_repeat_len ' + str(max_repeat_len) \
                                   + ' --thread ' + str(threads) \
                                   + ' --flanking_len ' + str(flanking_len) \
                                   + ' --tandem_region_cutoff ' + str(tandem_region_cutoff) \
                                   + ' --ref_index ' + str(ref_index) \
                                   + ' -r ' + reference
            os.system(coarse_boundary_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1.1: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
        resut_file = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1.2: determine fine-grained TIR')
            # ??????TIR?????????
            tir_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_TIR_transposons.py -g ' \
                                         + cut_reference + ' --seqs ' + longest_repeats_flanked_path\
                                         + ' -t ' + str(threads)+' --TRsearch_dir ' + TRsearch_dir \
                                         + ' --tmp_output_dir ' + tmp_output_dir \
                                         + ' --tandem_region_cutoff ' + str(tandem_region_cutoff) \
                                         + ' --ref_index ' + str(ref_index) \
                                         + ' --plant ' + str(plant) + ' --flanking_len ' + str(flanking_len)
            log.logger.debug(tir_identification_command)
            os.system(tir_identification_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1.2: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        resut_file = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1.3: determine fine-grained Helitron')
            # ??????Helitron?????????
            helitron_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Helitron_transposons.py --seqs ' \
                                              + longest_repeats_flanked_path + ' -g ' + cut_reference + ' -t ' + str(threads) \
                                              + ' --tmp_output_dir ' + tmp_output_dir + ' --EAHelitron ' + EAHelitron \
                                              + ' --ref_index ' + str(ref_index) + ' --flanking_len ' + str(flanking_len)

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
            log.logger.info("Running time of step1.3: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')


        resut_file = tmp_output_dir + '/confident_other_'+str(ref_index)+'.fa'
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1.4: homology-based other TE searching')
            #???????????????????????????
            other_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Other_transposons.py ' \
                                           + ' --seqs ' + longest_repeats_flanked_path\
                                           + ' -t ' + str(threads) \
                                           + ' --tmp_output_dir ' + tmp_output_dir + ' --query_coverage ' + str(0.8) \
                                           + ' --subject_coverage ' + str(0) + ' --ref_index ' + str(ref_index) \
                                           + ' --library_dir ' + str(library_dir)
            os.system(other_identification_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1.4: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')


    # ??????TIR??????????????????LTR????????????intact LTR or LTR terminals or LTR internals???
    # 1.1 ????????????parts???TIR??????
    confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    confident_helitron_path = tmp_output_dir + '/confident_helitron.fa'
    confident_other_path = tmp_output_dir + '/confident_other.fa'
    os.system('rm -f ' + confident_tir_path)
    os.system('rm -f ' + confident_helitron_path)
    os.system('rm -f ' + confident_other_path)
    for ref_index, ref_rename_path in enumerate(cut_references):
        cur_confident_tir_path = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
        cur_confident_helitron_path = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.fa'
        cur_confident_other_path = tmp_output_dir + '/confident_other_' + str(ref_index) + '.fa'
        os.system('cat ' + cur_confident_tir_path + ' >> ' + confident_tir_path)
        os.system('cat ' + cur_confident_helitron_path + ' >> ' + confident_helitron_path)
        os.system('cat ' + cur_confident_other_path + ' >> ' + confident_other_path)

    log.logger.info('Start step2: Structural Based LTR Searching')
    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    resut_file = confident_ltr_cut_path
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        LTR_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_LTR_transposons.py ' \
                                       + ' -g ' + reference + ' --ltrfinder_home ' + LTR_finder_parallel_Home \
                                       + ' -t ' + str(threads) \
                                       + ' --tmp_output_dir ' + tmp_output_dir \
                                       + ' --recover ' + str(recover) + ' --miu ' + str(miu)
        os.system(LTR_identification_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step2: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')


    # log.logger.info('Start step3: Remove nested TE')
    # starttime = time.time()
    # remove_nested_command = 'cd ' + test_home + ' && python3 ' + test_home + '/remove_nested.py ' \
    #                              + ' -g ' + reference + ' --confident_ltr_cut ' + confident_ltr_cut_path \
    #                              + ' --confident_tir ' + confident_tir_path \
    #                              + ' --confident_helitron ' + confident_helitron_path \
    #                              + ' --confident_other ' + confident_other_path \
    #                              + ' -t ' + str(threads) \
    #                              + ' --tmp_output_dir ' + tmp_output_dir \
    #                              + ' --global_flanking_filter ' + str(global_flanking_filter) \
    #                              + ' --remove_nested ' + str(remove_nested) + ' --test_home ' +str(test_home)

    # os.system(remove_nested_command)
    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.info("Running time of step3: %.8s s" % (dtime))

    # confident_TE_path = tmp_output_dir + '/confident_TE.fa'

    starttime = time.time()
    log.logger.info('Start step3: generate non-redundant library')
    generate_lib_command = 'cd ' + test_home + ' && python3 ' + test_home + '/get_nonRedundant_lib.py' \
                           + ' --confident_ltr_cut ' + confident_ltr_cut_path \
                           + ' --confident_tir ' + confident_tir_path \
                           + ' --confident_helitron ' + confident_helitron_path \
                           + ' --confident_other ' + confident_other_path \
                           + ' -t ' + str(threads) + ' --tmp_output_dir ' + tmp_output_dir

    os.system(generate_lib_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of step3: %.8s s" % (dtime))

    confident_TE_consensus = tmp_output_dir + '/confident_TE.cons.fa'
    starttime = time.time()
    log.logger.info('Start step4: generate classified library')
    TEClass_home = os.getcwd() + '/classification'
    classify_lib_command = 'cd ' + test_home + ' && python3 ' + test_home + '/get_classified_lib.py' \
                           + ' --confident_TE_consensus ' + confident_TE_consensus \
                           + ' -t ' + str(threads) + ' --tmp_output_dir ' + tmp_output_dir \
                           + ' --classified ' + str(classified) + ' --TEClass_home ' + str(TEClass_home) \
                           + ' --debug ' + str(debug)

    os.system(classify_lib_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of step4: %.8s s" % (dtime))

    starttime = time.time()
    log.logger.info('Start step5: clean library')
    clean_lib_command = 'cd ' + test_home + ' && python3 ' + test_home + '/clean_lib.py' \
                           + ' --tmp_output_dir ' + tmp_output_dir \
                           + ' --debug ' + str(debug)

    os.system(clean_lib_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of step5: %.8s s" % (dtime))


    pipeline_endtime = time.time()
    dtime = pipeline_endtime - pipeline_starttime
    log.logger.info("Running time of the whole pipeline: %.8s s" % (dtime))
