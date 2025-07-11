#-- coding: UTF-8 --
import argparse
import re
import shutil
import subprocess
import datetime
import json
import os
import sys
import time
from multiprocessing import cpu_count
from pathlib import Path
from module.Util import Logger, file_exist, read_fasta, filter_short_contigs_in_genome, create_or_clear_directory, \
    copy_files, update_prev_TE, get_fixed_extend_base_threshold, get_genome_size

current_folder = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.join(current_folder, ".")

if __name__ == '__main__':
    # We define the default parameters for HiTE.
    default_threads = int(cpu_count())
    default_fixed_extend_base_threshold = 4000
    default_chunk_size = 400
    default_tandem_region_cutoff = 0.5
    default_max_single_repeat_len = 30000
    default_plant = 1
    default_te_type = 'all'
    default_curated_lib = None
    default_recover = 0
    default_annotate = 0
    default_search_struct = 1
    default_BM_RM2 = 0
    default_BM_EDTA = 0
    default_BM_HiTE = 0
    default_EDTA_home = ''
    default_coverage_threshold = 0.95
    default_skip_HiTE = 0
    default_flanking_len = 50
    default_is_denovo_nonltr = 1
    default_debug = 0
    default_chrom_seg_length = 1_000_000
    default_classified = 1
    default_domain = 0
    default_miu = str(1.3e-8)
    default_remove_nested = 1
    default_use_HybridLTR = 1
    default_use_NeuralTE = 1
    default_is_wicker = 0
    default_is_output_LTR_lib = 1
    default_min_TE_len = 80

    version_num = '3.3.3'

    describe_image = '\n' + \
    '     __  __     __     ______   ______    \n' + \
    '    /\ \_\ \   /\ \   /\__  _\ /\  ___\   \n' + \
    '    \ \  __ \  \ \ \  \/_/\ \/ \ \  __\   \n' + \
    '     \ \_\ \_\  \ \_\    \ \_\  \ \_____\ \n' + \
    '      \/_/\/_/   \/_/     \/_/   \/_____/ ' + \
    'version ' + str(version_num) + '\n\n'
    print(describe_image)

    # 1.parse args
    describe_info = '########################## HiTE, version ' + str(version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--genome', required=True, metavar='genome', help='Input genome assembly path')
    parser.add_argument("--out_dir", nargs="?", default=os.getcwd(), help="The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.")

    parser.add_argument("--work_dir", nargs="?", default='/tmp', help="The temporary work directory for HiTE.")
    parser.add_argument('--thread', metavar='thread_num', help='Input thread num, default = [ '+str(default_threads)+' ]')
    parser.add_argument('--chunk_size', metavar='chunk_size', help='The chunk size of genome, default = [ ' + str(default_chunk_size) + ' MB ]')
    parser.add_argument('--miu', metavar='miu', help='The neutral mutation rate (per bp per ya), default = [ ' + str(default_miu) + ' ]')
    parser.add_argument('--plant', metavar='is_plant', help='Is it a plant genome, 1: true, 0: false. default = [ ' + str(default_plant) + ' ]')
    parser.add_argument('--te_type', metavar='te_type', help='Retrieve specific type of TE output [ltr|tir|helitron|non-ltr|all]. default = [ ' + str(default_te_type) + ' ]')
    parser.add_argument('--curated_lib', metavar='curated_lib', help='Provide a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome. We recommend using TE libraries from Repbase and ensuring the format follows >header#class_name. default = [ ' + str(default_curated_lib) + ' ]')
    # parser.add_argument('--classified', metavar='is_classified', help='Whether to classify TE models, HiTE uses RepeatClassifier from RepeatModeler to classify TEs, 1: true, 0: false. default = [ ' + str(default_classified) + ' ]')
    parser.add_argument('--remove_nested', metavar='is_remove_nested',help='Whether to remove nested TE, 1: true, 0: false. default = [ ' + str(default_remove_nested) + ' ]')
    parser.add_argument('--domain', metavar='is_domain', help='Whether to obtain TE domains, HiTE uses RepeatPeps.lib from RepeatMasker to obtain TE domains, 1: true, 0: false. default = [ ' + str(default_domain) + ' ]')
    parser.add_argument('--recover', metavar='is_recover', help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ ' + str(default_recover) + ' ]')
    parser.add_argument('--annotate', metavar='is_annotate', help='Whether to annotate the genome using the TE library generated, 1: true, 0: false. default = [ ' + str(default_annotate) + ' ]')
    parser.add_argument('--search_struct', metavar='search_struct', help='Is the structural information of full-length copies being searched, 1: true, 0: false. default = [ ' + str(default_search_struct) + ' ]')
    parser.add_argument('--BM_RM2', metavar='BM_RM2', help='Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false. default = [ ' + str(default_BM_RM2) + ' ]')
    parser.add_argument('--BM_EDTA', metavar='BM_EDTA', help='Whether to conduct benchmarking of EDTA, 1: true, 0: false. default = [ ' + str(default_BM_EDTA) + ' ]')
    parser.add_argument('--BM_HiTE', metavar='BM_HiTE', help='Whether to conduct benchmarking of HiTE, 1: true, 0: false. default = [ ' + str(default_BM_HiTE) + ' ]')
    parser.add_argument('--EDTA_home', metavar='EDTA_home', help='When conducting benchmarking of EDTA, you will be asked to input EDTA home path.')
    parser.add_argument('--coverage_threshold', metavar='coverage_threshold', help='The coverage threshold of benchmarking methods.')
    parser.add_argument('--species', metavar='species', help='Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).')
    parser.add_argument('--skip_HiTE', metavar='skip_HiTE', help='Whether to skip_HiTE, 1: true, 0: false. default = [ ' + str(default_skip_HiTE) + ' ]')
    parser.add_argument('--is_denovo_nonltr', metavar='is_denovo_nonltr', help='Whether to detect non-ltr de novo, 1: true, 0: false. default = [ ' + str(default_is_denovo_nonltr) + ' ]')
    parser.add_argument('--debug', metavar='is_debug', help='Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ ' + str(default_debug) + ' ]')
    parser.add_argument('--use_HybridLTR', metavar='use_HybridLTR', help='Whether to use HybridLTR to identify LTRs, 1: true, 0: false. default = [' + str(default_use_HybridLTR) + ' ]')
    parser.add_argument('--use_NeuralTE', metavar='use_NeuralTE', help='Whether to use NeuralTE to classify TEs, 1: true, 0: false. default = [' + str(default_use_NeuralTE) + ' ]')
    parser.add_argument('--is_wicker', metavar='is_wicker', help='Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [ ' + str(default_is_wicker) + ' ]')
    parser.add_argument('--is_output_LTR_lib', metavar='is_output_LTR_lib', help='Whether to output LTR library. default = [ ' + str(default_is_output_LTR_lib) + ' ]')
    parser.add_argument('--min_TE_len', metavar='min_TE_len', help='The minimum TE length, default = [ ' + str(default_min_TE_len) + ' bp ]')

    parser.add_argument('--flanking_len', metavar='flanking_len', help='The flanking length of candidates to find the true boundaries, default = [ ' + str(default_flanking_len) + ' ]')
    parser.add_argument('--fixed_extend_base_threshold', metavar='fixed_extend_base_threshold', help='The length of variation can be tolerated during pairwise alignment, default = [ '+str(default_fixed_extend_base_threshold)+' ]')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff', help='Cutoff of the candidates regarded as tandem region, default = [ '+str(default_tandem_region_cutoff)+' ]')
    parser.add_argument('--max_repeat_len', metavar='max_repeat_len', help='The maximum length of a single repeat, default = [ ' + str(default_max_single_repeat_len) + ' ]')
    parser.add_argument('--chrom_seg_length', metavar='chrom_seg_length', help='The length of genome segments, default = [ ' + str(default_chrom_seg_length) + ' ]')
    parser.add_argument('--shared_prev_TE', metavar='shared_prev_TE', help='The path of shared previous TEs')

    args = parser.parse_args()

    reference = args.genome
    threads = args.thread
    fixed_extend_base_threshold = args.fixed_extend_base_threshold
    chunk_size = args.chunk_size
    tandem_region_cutoff = args.tandem_region_cutoff
    max_repeat_len = args.max_repeat_len
    chrom_seg_length = args.chrom_seg_length
    flanking_len = args.flanking_len
    plant = args.plant
    te_type = args.te_type
    curated_lib = args.curated_lib
    # classified = args.classified
    remove_nested = args.remove_nested
    domain = args.domain
    miu = args.miu
    recover = args.recover
    annotate = args.annotate
    search_struct = args.search_struct
    BM_RM2 = args.BM_RM2
    BM_EDTA = args.BM_EDTA
    BM_HiTE = args.BM_HiTE
    EDTA_home = args.EDTA_home
    coverage_threshold = args.coverage_threshold
    species = args.species
    skip_HiTE = args.skip_HiTE
    is_denovo_nonltr = args.is_denovo_nonltr
    debug = args.debug
    use_HybridLTR = args.use_HybridLTR
    use_NeuralTE = args.use_NeuralTE
    is_wicker = args.is_wicker
    is_output_LTR_lib = args.is_output_LTR_lib
    min_TE_len = args.min_TE_len
    shared_prev_TE = args.shared_prev_TE

    output_dir = os.path.abspath(args.out_dir)
    os.makedirs(output_dir, exist_ok=True)
    work_dir = os.path.abspath(args.work_dir)
    Path(work_dir).mkdir(exist_ok=True)

    if reference is None:
        print('\nreference path can not be empty')
        parser.print_help()
        exit(-1)

    if not os.path.exists(reference):
        print('\nCannot find input genome assembly: ' + str(reference))
        parser.print_help()
        exit(-1)

    if not os.path.isabs(reference):
        reference = os.path.abspath(reference)

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    if use_HybridLTR is None:
        use_HybridLTR = int(default_use_HybridLTR)
    else:
        use_HybridLTR = int(use_HybridLTR)

    if use_NeuralTE is None:
        use_NeuralTE = int(default_use_NeuralTE)
    else:
        use_NeuralTE = int(use_NeuralTE)

    if is_wicker is None:
        is_wicker = int(default_is_wicker)
    else:
        is_wicker = int(is_wicker)

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

    if te_type is None:
        te_type = default_te_type
    else:
        te_type = te_type

    if curated_lib is None:
        curated_lib = default_curated_lib
    else:
        curated_lib = curated_lib

    if is_output_LTR_lib is None:
        is_output_LTR_lib = default_is_output_LTR_lib
    else:
        is_output_LTR_lib = int(is_output_LTR_lib)

    if min_TE_len is None:
        min_TE_len = default_min_TE_len
    else:
        min_TE_len = int(min_TE_len)

    if remove_nested is None:
        remove_nested = default_remove_nested
    else:
        remove_nested = int(remove_nested)

    if domain is None:
        domain = default_domain
    else:
        domain = int(domain)
    
    if miu is None:
        miu = default_miu
    else:
        miu = str(miu)

    if recover is None:
        recover = default_recover
    else:
        recover = int(recover)

    if annotate is None:
        annotate = default_annotate
    else:
        annotate = int(annotate)

    if search_struct is None:
        search_struct = default_search_struct
    else:
        search_struct = int(search_struct)

    if BM_RM2 is None:
        BM_RM2 = default_BM_RM2
    else:
        BM_RM2 = int(BM_RM2)

    if BM_EDTA is None:
        BM_EDTA = default_BM_EDTA
    else:
        BM_EDTA = int(BM_EDTA)

    if BM_HiTE is None:
        BM_HiTE = default_BM_HiTE
    else:
        BM_HiTE = int(BM_HiTE)

    if EDTA_home is None:
        EDTA_home = default_EDTA_home
    else:
        EDTA_home = str(EDTA_home)

    if coverage_threshold is None:
        coverage_threshold = default_coverage_threshold
    else:
        coverage_threshold = float(coverage_threshold)

    if skip_HiTE is None:
        skip_HiTE = default_skip_HiTE
    else:
        skip_HiTE = int(skip_HiTE)

    if is_denovo_nonltr is None:
        is_denovo_nonltr = default_is_denovo_nonltr
    else:
        is_denovo_nonltr = int(is_denovo_nonltr)

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

    total_starttime = time.time()
    tools_dir = project_dir + '/tools'

    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

    # 创建本地临时目录，存储计算结果
    tmp_output_dir = output_dir
    if recover == 0:
        # 如果目录已存在，则删除该目录及其中的所有内容
        if os.path.exists(tmp_output_dir):
            shutil.rmtree(tmp_output_dir)
    # 创建新的目录
    os.makedirs(tmp_output_dir, exist_ok=True)

    log = Logger(tmp_output_dir + '/HiTE.log', level='debug')

    try:
        shutil.copy2(reference, os.path.join(tmp_output_dir, os.path.basename(reference)))
    except Exception as e:
        log.logger.error(f"Error copying file: {e}")

    raw_reference = os.path.join(tmp_output_dir, os.path.basename(reference))
    genome_size = get_genome_size(raw_reference)
    fixed_extend_base_threshold = get_fixed_extend_base_threshold(genome_size)

    all_low_copy_tir = os.path.join(tmp_output_dir, 'tir_low_copy.fa')
    all_low_copy_helitron = os.path.join(tmp_output_dir, 'helitron_low_copy.fa')
    all_low_copy_non_ltr = os.path.join(tmp_output_dir, 'non_ltr_low_copy.fa')
    if os.path.exists(all_low_copy_tir):
        os.remove(all_low_copy_tir)
    if os.path.exists(all_low_copy_helitron):
        os.remove(all_low_copy_helitron)
    if os.path.exists(all_low_copy_non_ltr):
        os.remove(all_low_copy_non_ltr)

    # preset steps, no needs to execute time-consuming job again
    # when the job is retried due to abnormal termination.
    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    all_te_types = ['ltr', 'tir', 'helitron', 'non-ltr', 'all']
    if te_type not in all_te_types:
        log.logger.error('Specified an invalid TE type: ' + te_type + '. Please choose from ' + str(all_te_types))
        sys.exit(-1)

    if curated_lib is not None:
        curated_lib = os.path.realpath(curated_lib)
        curated_names, curated_contigs = read_fasta(curated_lib)
        if len(curated_names) <= 0:
            log.logger.error('You have provided an invalid or empty curated library: ' + str(curated_lib) + '. Please check and correct it.')
            sys.exit(-1)

    if BM_EDTA == 1 and not os.path.exists(EDTA_home + '/lib-test.pl'):
        print('Cannot conduct benchmarking of EDTA, Invalid EDTA home: ' + EDTA_home)
        sys.exit(-1)

    log.logger.info('\n-------------------------------------------------------------------------------------------\n'
                    'Copyright (C) 2022 Kang Hu ( kanghu@csu.edu.cn )\n'
                    'Hunan Provincial Key Lab on Bioinformatics, School of Computer Science and \n'
                    'Engineering, Central South University, Changsha 410083, P.R. China.\n'
                    '-------------------------------------------------------------------------------------------')

    log.logger.info('\nParameters configuration\n'
                    '====================================System settings========================================\n'
                    '  [Setting] Reference sequences / assemblies path = [ ' + str(reference) + ' ]\n'
                    '  [Setting] Temporary work dir = [ ' + str(work_dir) + ' ] Default( ' + str('/tmp') + ' )\n'
                    '  [Setting] Is remove nested TE = [ ' + str(remove_nested) + ' ] Default( ' + str(default_remove_nested) + ' )\n'
                    '  [Setting] Is getting domain = [ ' + str(domain) + ' ] Default( ' + str(default_domain) + ' )\n'
                    '  [Setting] The neutral mutation rate (per bp per ya) = [ ' + str(miu) + ' ] Default( ' + str(default_miu) + ' )\n'
                    '  [Setting] Threads = [ ' + str(threads) + ' ]  Default( ' + str(default_threads) + ' )\n'
                    '  [Setting] The chunk size of large genome = [ ' + str(chunk_size) + ' ] MB Default( ' + str(default_chunk_size) + ' ) MB\n'
                    '  [Setting] Is plant genome = [ ' + str(plant) + ' ]  Default( ' + str(default_plant) + ' )\n'
                    '  [Setting] Retrieve specific type of TE output = [ ' + str(te_type) + ' ]  Default( ' + str(default_te_type) + ' )\n'
                    '  [Setting] Curated library = [ ' + str(curated_lib) + ' ]  Default( ' + str(default_curated_lib) + ' )\n'
                    '  [Setting] recover = [ ' + str(recover) + ' ]  Default( ' + str(default_recover) + ' )\n'
                    '  [Setting] annotate = [ ' + str(annotate) + ' ]  Default( ' + str(default_annotate) + ' )\n'                                                                                    
                    '  [Setting] search_struct = [ ' + str(search_struct) + ' ] Default( ' + str(default_search_struct) + ' )\n'                                                                                        
                    '  [Setting] BM_RM2 = [ ' + str(BM_RM2) + ' ]  Default( ' + str(default_BM_RM2) + ' )\n'
                    '  [Setting] BM_EDTA = [ ' + str(BM_EDTA) + ' ]  Default( ' + str(default_BM_EDTA) + ' )\n'
                    '  [Setting] BM_HiTE = [ ' + str(BM_HiTE) + ' ]  Default( ' + str(default_BM_HiTE) + ' )\n'
                    '  [Setting] EDTA_home = [' + str(EDTA_home) + ']\n'
                    '  [Setting] coverage_threshold = [ ' + str(coverage_threshold) + ' ]  Default( ' + str(default_coverage_threshold) + ' )\n'
                    '  [Setting] skip_HiTE = [ ' + str(skip_HiTE) + ' ]  Default( ' + str(default_skip_HiTE) + ' )\n'
                    '  [Setting] is_denovo_nonltr = [ ' + str(is_denovo_nonltr) + ' ]  Default( ' + str(default_is_denovo_nonltr) + ' )\n'
                    '  [Setting] use_HybridLTR = [ ' + str(use_HybridLTR) + ' ]  Default( ' + str(default_use_HybridLTR) + ' )\n'                                                                                                                                    
                    '  [Setting] use_NeuralTE = [ ' + str(use_NeuralTE) + ' ]  Default( ' + str(default_use_NeuralTE) + ' )\n'
                    '  [Setting] is_wicker = [ ' + str(is_wicker) + ' ]  Default( ' + str(default_is_wicker) + ' )\n'
                    '  [Setting] debug = [ ' + str(debug) + ' ]  Default( ' + str(default_debug) + ' )\n'
                    '  [Setting] Output Directory = [' + str(output_dir) + ']\n'
                                                                                                                                                                                                           
                    '  [Setting] Fixed extend bases threshold = [ ' + str(fixed_extend_base_threshold) + ' ] Default( ' + str(default_fixed_extend_base_threshold) + ' )\n'
                    '  [Setting] Flanking length of TE = [ ' + str(flanking_len) + ' ]  Default( ' + str(default_flanking_len) + ' )\n'
                    '  [Setting] Cutoff of the repeat regarded as tandem sequence = [ ' + str(tandem_region_cutoff) + ' ] Default( ' + str(default_tandem_region_cutoff) + ' )\n'
                    '  [Setting] The length of genome segments = [ ' + str(chrom_seg_length) + ' ]  Default( ' + str(default_chrom_seg_length) + ' )'
                    )

    TRsearch_dir = tools_dir
    test_home = project_dir + '/module'

    # The HiTE pipeline performs recognition of LTR, Non-LTR, TIR, and Helitron transposons.
    # The organizational structure of HiTE is as follows:
    # Pipeline: main.py
    #     ├──LTR: judge_LTR_transposons.py
    #     ├──Homology-Non-LTR: judge_Other_transposons.py
    #     ├──split genome into chunks: split_genome_chunks.py
    #       ├──De novo TE searching: coarse_boundary.py
    #       ├──TIR: judge_TIR_transposons.py
    #       ├──Helitron: judge_Helitron_transposons.py
    #       └──De novo-Non-LTR: judge_Non_LTR_transposons.py
    #     ├──generate TE library: get_nonRedundant_lib.py
    #       └──unwrap nested TE: remove_nested_lib.py
    #     ├──genome annotation: annotate_genome.py
    #     ├──benchmarking reproduction: benchmarking.py
    #     └──clean temporary files: clean_lib.py

    if threads < 4:
        log.logger.warning('The number of threads running HiTE is less than 4, we strongly recommend using a higher number of threads.')
    # reduce the threads number to avoid run out of resources
    threads = max(1, threads-4)

    pipeline_starttime = time.time()
    if skip_HiTE != 1:
        reference_clean = tmp_output_dir + '/genome.fa.clean'
        resut_file = reference_clean
        if not is_recover or not file_exist(resut_file):
            starttime = time.time()
            log.logger.info('Start step1: Remove redundant contigs from a genome assembly using minimap2')
            genome_clean_command = 'genome_clean.py ' \
                                           + ' -i ' + raw_reference \
                                           + ' -o ' + reference_clean \
                                           + ' -t ' + str(threads) \
                                           + ' -w ' + str(work_dir)
            log.logger.info(genome_clean_command)
            os.system(genome_clean_command)
            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of step1: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')
        reference = reference_clean

        confident_other_path = tmp_output_dir + '/confident_other.fa'
        resut_file = confident_other_path
        if not is_recover or not file_exist(resut_file):
            if te_type == 'all' or te_type == 'non-ltr':
                starttime = time.time()
                log.logger.info('Start step2: homology-based other TE searching')
                other_identification_command = 'judge_Other_transposons.py ' \
                                               + ' -r ' + reference \
                                               + ' -t ' + str(threads) \
                                               + ' --tmp_output_dir ' + tmp_output_dir  \
                                               + ' --recover ' + str(recover) \
                                               + ' --min_TE_len ' + str(min_TE_len) \
                                               + ' -w ' + str(work_dir)
                log.logger.info(other_identification_command)
                os.system(other_identification_command)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step2: %.8s s" % (dtime))
            else:
                # 创建一个空的输出文件，以跳过nextflow的检查
                confident_other = os.path.join(tmp_output_dir, "confident_other.fa")
                empty_files = [confident_other]
                for empty_file in empty_files:
                    os.system('touch ' + empty_file)
        else:
            log.logger.info(resut_file + ' exists, skip...')

        # --------------------------------------------------------------------------------------
        starttime = time.time()
        log.logger.info('Start step3.0: Splitting genome assembly into chunks')
        split_genome_command = 'split_genome_chunks.py -g ' \
                                     + reference + ' --tmp_output_dir ' + tmp_output_dir \
                                     + ' --chrom_seg_length ' + str(chrom_seg_length) + ' --chunk_size ' + str(chunk_size)
        log.logger.info(split_genome_command)
        os.system(split_genome_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step3.0: %.8s s" % (dtime))

        reg_str = 'genome.cut(\d+).fa$'
        cut_references = []
        for filename in os.listdir(tmp_output_dir):
            match = re.search(reg_str, filename)
            if match:
                ref_index = match.group(1)
                cut_references.append((ref_index, tmp_output_dir + '/' + filename))

        # Using identified TEs to mask the genome in order to reduce computational load in all-vs-all alignments.
        if shared_prev_TE is not None:
            prev_TE = shared_prev_TE
        else:
            prev_TE = tmp_output_dir + '/prev_TE.fa'
        if curated_lib is not None and os.path.exists(curated_lib):
            update_prev_TE(prev_TE, curated_lib)
        # The outcomes of homologous methods can only serve as supplementary information and should not be used as masks,
        # as this could potentially obscure many genuine non-LTR local masks, rendering the de novo method unable to identify them.
        # os.system('cat ' + confident_other_path + ' >> ' + prev_TE)

        split_ref_dir = tmp_output_dir + '/ref_chr'
        for cut_reference_item in cut_references:
            ref_index = cut_reference_item[0]
            cut_reference = cut_reference_item[1]
            log.logger.info('Current chunk: ' + str(ref_index))

            longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
            longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
            resut_file = longest_repeats_path
            if not is_recover or not file_exist(resut_file) or not file_exist(longest_repeats_flanked_path):
                if te_type != 'ltr':
                    starttime = time.time()
                    log.logger.info('Start 3.1: Coarse-grained boundary mapping')
                    coarse_boundary_command = 'coarse_boundary.py ' \
                                           + ' -g ' + cut_reference + ' --tmp_output_dir ' + tmp_output_dir \
                                           + ' --prev_TE ' + str(prev_TE) \
                                           + ' --fixed_extend_base_threshold ' + str(fixed_extend_base_threshold) \
                                           + ' --max_repeat_len ' + str(max_repeat_len) \
                                           + ' --thread ' + str(threads) \
                                           + ' --flanking_len ' + str(flanking_len) \
                                           + ' --tandem_region_cutoff ' + str(tandem_region_cutoff) \
                                           + ' --ref_index ' + str(ref_index) \
                                           + ' -r ' + reference + ' --recover ' + str(recover) \
                                           + ' --debug ' + str(debug) + ' -w ' + str(work_dir)
                    log.logger.info(coarse_boundary_command)
                    os.system(coarse_boundary_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step3.1: %.8s s" % (dtime))
            else:
                log.logger.info(resut_file + ' exists, skip...')

            longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
            resut_file = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
            if not is_recover or not file_exist(resut_file):
                if te_type == 'all' or te_type == 'tir':
                    starttime = time.time()
                    log.logger.info('Start step3.2: determine fine-grained TIR')
                    tir_identification_command = 'judge_TIR_transposons.py ' \
                                                 + ' --seqs ' + longest_repeats_flanked_path \
                                                 + ' -t ' + str(threads) \
                                                 + ' --tmp_output_dir ' + tmp_output_dir \
                                                 + ' --tandem_region_cutoff ' + str(tandem_region_cutoff) \
                                                 + ' --ref_index ' + str(ref_index) \
                                                 + ' --plant ' + str(plant) \
                                                 + ' --flanking_len ' + str(flanking_len) \
                                                 + ' --recover ' + str(recover) \
                                                 + ' --debug ' + str(debug) \
                                                 + ' -r ' + reference \
                                                 + ' --split_ref_dir ' + split_ref_dir \
                                                 + ' --prev_TE  ' + prev_TE \
                                                 + ' --all_low_copy_tir ' + all_low_copy_tir \
                                                 + ' --min_TE_len ' + str(min_TE_len) \
                                                 + ' -w ' + str(work_dir)
                    log.logger.debug(tir_identification_command)
                    os.system(tir_identification_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step3.2: %.8s s" % (dtime))
                else:
                    # 创建一个空的输出文件，以跳过nextflow的检查
                    confident_tir = os.path.join(tmp_output_dir, 'confident_tir_'+str(ref_index)+'.fa')
                    empty_files = [confident_tir]
                    for empty_file in empty_files:
                        os.system('touch ' + empty_file)
            else:
                log.logger.info(resut_file + ' exists, skip...')

            resut_file = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
            if not is_recover or not file_exist(resut_file):
                if te_type == 'all' or te_type == 'helitron':
                    starttime = time.time()
                    log.logger.info('Start step3.3: determine fine-grained Helitron')
                    helitron_identification_command = 'judge_Helitron_transposons.py --seqs ' \
                                                      + longest_repeats_flanked_path + ' -r ' + reference + ' -t ' + str(threads) \
                                                      + ' --tmp_output_dir ' + tmp_output_dir \
                                                      + ' --ref_index ' + str(ref_index) \
                                                      + ' --flanking_len ' + str(flanking_len) \
                                                      + ' --recover ' + str(recover) \
                                                      + ' --debug ' + str(debug) \
                                                      + ' --split_ref_dir ' + split_ref_dir \
                                                      + ' --prev_TE  ' + prev_TE \
                                                      + ' --all_low_copy_helitron ' + all_low_copy_helitron \
                                                      + ' --min_TE_len ' + str(min_TE_len) \
                                                      + ' -w ' + str(work_dir)

                    log.logger.info(helitron_identification_command)
                    os.system(helitron_identification_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step3.3: %.8s s" % (dtime))
                else:
                    # 创建一个空的输出文件，以跳过nextflow的检查
                    confident_helitron = os.path.join(tmp_output_dir, 'confident_helitron_' + str(ref_index) + '.fa')
                    empty_files = [confident_helitron]
                    for empty_file in empty_files:
                        os.system('touch ' + empty_file)
            else:
                log.logger.info(resut_file + ' exists, skip...')

            resut_file = tmp_output_dir + '/confident_non_ltr_' + str(ref_index) + '.fa'
            if not is_recover or not file_exist(resut_file):
                if te_type == 'all' or te_type == 'non-ltr':
                    starttime = time.time()
                    log.logger.info('Start step3.4: determine fine-grained Non-LTR')
                    non_ltr_identification_command = 'judge_Non_LTR_transposons.py'\
                                                 + ' --seqs ' + longest_repeats_flanked_path + ' -t ' + str(threads) \
                                                 + ' --tmp_output_dir ' + tmp_output_dir \
                                                 + ' --recover ' + str(recover) \
                                                 + ' --plant ' + str(plant) \
                                                 + ' --debug ' + str(debug) \
                                                 + ' --flanking_len ' + str(flanking_len) \
                                                 + ' --ref_index ' + str(ref_index) \
                                                 + ' --split_ref_dir ' + split_ref_dir \
                                                 + ' --is_denovo_nonltr ' + str(is_denovo_nonltr) \
                                                 + ' -r ' + reference \
                                                 + ' --prev_TE  ' + prev_TE \
                                                 + ' --all_low_copy_non_ltr ' + all_low_copy_non_ltr \
                                                 + ' --min_TE_len ' + str(min_TE_len) \
                                                 + ' -w ' + str(work_dir)
                    log.logger.debug(non_ltr_identification_command)
                    os.system(non_ltr_identification_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step3.4: %.8s s" % (dtime))
                else:
                    # 创建一个空的输出文件，以跳过nextflow的检查
                    confident_non_ltr = os.path.join(tmp_output_dir, 'confident_non_ltr_' + str(ref_index) + '.fa')
                    empty_files = [confident_non_ltr]
                    for empty_file in empty_files:
                        os.system('touch ' + empty_file)
            else:
                log.logger.info(resut_file + ' exists, skip...')

        # Merging all TEs generated from individual chunks.
        confident_tir_path = tmp_output_dir + '/confident_tir_merge.fa'
        confident_helitron_path = tmp_output_dir + '/confident_helitron_merge.fa'
        confident_non_ltr_path = tmp_output_dir + '/confident_non_ltr_merge.fa'
        tmp_lib = tmp_output_dir + '/tmp_lib.fa'
        if os.path.exists(confident_tir_path):
            os.remove(confident_tir_path)
        if os.path.exists(confident_helitron_path):
            os.remove(confident_helitron_path)
        if os.path.exists(confident_non_ltr_path):
            os.remove(confident_non_ltr_path)
        for ref_index, ref_rename_path in enumerate(cut_references):
            cur_confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.fa'
            cur_confident_helitron_path = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.fa'
            cur_confident_non_ltr_path = tmp_output_dir + '/confident_non_ltr_' + str(ref_index) + '.fa'
            if os.path.exists(cur_confident_tir_path):
                os.system('cat ' + cur_confident_tir_path + ' >> ' + confident_tir_path)
            if os.path.exists(cur_confident_helitron_path):
                os.system('cat ' + cur_confident_helitron_path + ' >> ' + confident_helitron_path)
            if os.path.exists(cur_confident_non_ltr_path):
                os.system('cat ' + cur_confident_non_ltr_path + ' >> ' + confident_non_ltr_path)
        if curated_lib is not None and os.path.exists(curated_lib):
            os.system('cat ' + curated_lib + ' > ' + tmp_lib)
        os.system('cat ' + confident_tir_path + ' ' + confident_helitron_path + ' ' + confident_non_ltr_path + ' ' + confident_other_path + ' >> ' + tmp_lib)

        confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
        resut_file = confident_ltr_cut_path
        if not is_recover or not file_exist(resut_file):
            if te_type == 'all' or te_type == 'ltr':
                log.logger.info('Start step4: Structural Based LTR Searching')
                starttime = time.time()
                TEClass_home = project_dir + '/classification'
                LTR_identification_command = 'judge_LTR_transposons.py ' \
                                             + ' -g ' + reference + ' -t ' + str(threads) \
                                             + ' --tmp_output_dir ' + tmp_output_dir + ' --recover ' + str(recover) \
                                             + ' --miu ' + str(miu) + ' --use_HybridLTR ' + str(use_HybridLTR) \
                                             + ' --use_NeuralTE ' + str(use_NeuralTE) + ' --is_wicker ' + str(is_wicker) \
                                             + ' --is_output_lib ' + str(is_output_LTR_lib) + ' --debug ' + str(debug) \
                                             + ' -w ' + str(work_dir) + ' --prev_TE ' + tmp_lib
                log.logger.info(LTR_identification_command)
                os.system(LTR_identification_command)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step4: %.8s s" % (dtime))
            else:
                # 创建一个空的输出文件，以跳过nextflow的检查
                confident_ltr_cut = os.path.join(tmp_output_dir, "confident_ltr_cut.fa")
                confident_ltr_internal = os.path.join(tmp_output_dir, "confident_ltr.internal.fa")
                confident_ltr_terminal = os.path.join(tmp_output_dir, "confident_ltr.terminal.fa")
                intact_LTR = os.path.join(tmp_output_dir, "intact_LTR.fa")
                classified_intact_LTR = os.path.join(tmp_output_dir, "intact_LTR.fa.classified")
                intact_LTR_list = os.path.join(tmp_output_dir, "intact_LTR.list")
                empty_files = [confident_ltr_cut, confident_ltr_internal, confident_ltr_terminal,
                               intact_LTR, classified_intact_LTR, intact_LTR_list]
                for empty_file in empty_files:
                    os.system('touch ' + empty_file)
        else:
            log.logger.info(resut_file + ' exists, skip...')

        # # tir, helitron, non_ltr detected by LTR module
        # tir_from_ltr_path = os.path.join(tmp_output_dir, "confident_tir_from_ltr.fa")
        # helitron_from_ltr_path = os.path.join(tmp_output_dir, "confident_helitron_from_ltr.fa")
        # # non_ltr_from_ltr_path = os.path.join(tmp_output_dir, "confident_non_ltr_from_ltr.fa")
        # if os.path.exists(tir_from_ltr_path):
        #     os.system('cat ' + tir_from_ltr_path + ' >> ' + confident_tir_path)
        # if os.path.exists(helitron_from_ltr_path):
        #     os.system('cat ' + helitron_from_ltr_path + ' >> ' + confident_helitron_path)

        starttime = time.time()
        log.logger.info('Start step5: generate non-redundant library')
        TEClass_home = project_dir + '/classification'
        generate_lib_command = 'get_nonRedundant_lib.py' \
                               + ' --confident_ltr_cut ' + confident_ltr_cut_path \
                               + ' --confident_tir ' + confident_tir_path \
                               + ' --confident_helitron ' + confident_helitron_path \
                               + ' --confident_non_ltr ' + confident_non_ltr_path \
                               + ' --confident_other ' + confident_other_path \
                               + ' -t ' + str(threads) + ' --tmp_output_dir ' + tmp_output_dir \
                               + ' --use_NeuralTE ' + str(use_NeuralTE) \
                               + ' --is_wicker ' + str(is_wicker) \
                               + ' --domain ' + str(domain) \
                               + ' --curated_lib ' + str(curated_lib) \
                               + ' -w ' + str(work_dir)
        log.logger.info(generate_lib_command)
        os.system(generate_lib_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step5: %.8s s" % (dtime))

    # merge low copy TEs
    low_confident_TEs_merge = os.path.join(tmp_output_dir, 'low_confident_TE.fa')
    if os.path.exists(low_confident_TEs_merge):
        os.remove(low_confident_TEs_merge)
    low_confident_TEs = os.path.join(tmp_output_dir, 'low_confident_TE.cons.fa')
    if os.path.exists(all_low_copy_tir):
        os.system('cat ' + all_low_copy_tir + ' >> ' + low_confident_TEs_merge)
    if os.path.exists(all_low_copy_helitron):
        os.system('cat ' + all_low_copy_helitron + ' >> ' + low_confident_TEs_merge)
    if os.path.exists(all_low_copy_non_ltr):
        os.system('cat ' + all_low_copy_non_ltr + ' >> ' + low_confident_TEs_merge)
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -d 0 -G 0 -g 1 -A 80 -i ' + low_confident_TEs_merge + ' -o ' + low_confident_TEs + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')

    # merge low confident and confident TEs
    confident_TE_consensus = os.path.join(tmp_output_dir, 'confident_TE.cons.fa')
    all_TEs = os.path.join(tmp_output_dir, 'all_TE.fa')
    if os.path.exists(all_TEs):
        os.remove(all_TEs)
    if os.path.exists(confident_TE_consensus):
        os.system('cat ' + confident_TE_consensus + ' >> ' + all_TEs)
    if os.path.exists(low_confident_TEs):
        os.system('cat ' + low_confident_TEs + ' >> ' + all_TEs)


    if not file_exist(confident_TE_consensus):
        log.logger.error('Error, Cannot find TE path: ' + confident_TE_consensus)
    else:
        starttime = time.time()
        log.logger.info('Start step6: annotate genome')
        TEClass_home = project_dir + '/classification'
        annotate_genome_command = 'pan_annotate_genome.py --threads ' + str(threads) \
                                  + ' --panTE_lib ' + confident_TE_consensus + ' --genome_name HiTE ' \
                                  + ' --reference ' + raw_reference + ' --annotate ' + str(annotate) \
                                  + ' --output_dir ' + tmp_output_dir + ' -w ' + str(work_dir)
        log.logger.info(annotate_genome_command)
        os.system(annotate_genome_command)

        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step6: %.8s s" % (dtime))

        starttime = time.time()
        log.logger.info('Start step7: Start conduct benchmarking of RepeatModeler2, EDTA, and HiTE')
        benchmarking_command = 'benchmarking.py' \
                            + ' --tmp_output_dir ' + tmp_output_dir \
                            + ' --BM_RM2 ' + str(BM_RM2) + ' --BM_EDTA ' + str(BM_EDTA) + ' --BM_HiTE ' + str(BM_HiTE) \
                            + ' --coverage_threshold ' + str(coverage_threshold) + ' -t ' + str(threads) \
                            + ' --TE_lib ' + str(confident_TE_consensus) \
                            + ' -r ' + raw_reference + ' --recover ' + str(recover) + ' -w ' + str(work_dir)
        if EDTA_home is not None and EDTA_home.strip() != '':
            benchmarking_command += ' --EDTA_home ' + str(EDTA_home)
        if species is None or species.strip() == '':
            benchmarking_command += ' --species test'
        else:
            benchmarking_command += ' --species ' + str(species)
        log.logger.info(benchmarking_command)
        os.system(benchmarking_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step7: %.8s s" % (dtime))


    pipeline_endtime = time.time()
    dtime = pipeline_endtime - pipeline_starttime
    log.logger.info("Running time of the whole pipeline: %.8s s" % (dtime))


    clean_lib_command = 'clean_lib.py' \
                           + ' --tmp_output_dir ' + tmp_output_dir \
                           + ' --debug ' + str(debug)
    os.system(clean_lib_command)