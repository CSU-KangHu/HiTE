#-- coding: UTF-8 --
import argparse
import re
import subprocess
import datetime
import json
import os
import sys
import time
from multiprocessing import cpu_count

from module.Util import Logger, file_exist, read_fasta

current_folder = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.join(current_folder, ".")

if __name__ == '__main__':
    # We define the default parameters for HiTE.
    default_threads = int(cpu_count())
    default_fixed_extend_base_threshold = 1000
    default_chunk_size = 400
    default_tandem_region_cutoff = 0.5
    default_max_single_repeat_len = 30000
    default_plant = 1
    default_te_type = 'all'
    default_curated_lib = None
    default_recover = 0
    default_annotate = 0
    default_intact_annotate = 0
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
    default_chrom_seg_length = 100000
    default_classified = 1
    default_domain = 0
    default_miu = str(1.3e-8)
    default_remove_nested = 1
    default_use_NeuralTE = 1
    default_is_wicker = 0

    version_num = '3.2'

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
    parser.add_argument('--outdir', required=True, metavar='output_dir', help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')

    parser.add_argument('--thread', metavar='thread_num', help='Input thread num, default = [ '+str(default_threads)+' ]')
    parser.add_argument('--chunk_size', metavar='chunk_size', help='The chunk size of genome, default = [ ' + str(default_chunk_size) + ' MB ]')
    parser.add_argument('--miu', metavar='miu', help='The neutral mutation rate (per bp per ya), default = [ ' + str(default_miu) + ' ]')
    parser.add_argument('--plant', metavar='is_plant', help='Is it a plant genome, 1: true, 0: false. default = [ ' + str(default_plant) + ' ]')
    parser.add_argument('--te_type', metavar='te_type', help='Retrieve specific type of TE output [ltr|tir|helitron|non-ltr|all]. default = [ ' + str(default_te_type) + ' ]')
    parser.add_argument('--curated_lib', metavar='curated_lib', help='Provide a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome. We recommend using TE libraries from Repbase. default = [ ' + str(default_curated_lib) + ' ]')
    # parser.add_argument('--classified', metavar='is_classified', help='Whether to classify TE models, HiTE uses RepeatClassifier from RepeatModeler to classify TEs, 1: true, 0: false. default = [ ' + str(default_classified) + ' ]')
    parser.add_argument('--remove_nested', metavar='is_remove_nested',help='Whether to remove nested TE, 1: true, 0: false. default = [ ' + str(default_remove_nested) + ' ]')
    parser.add_argument('--domain', metavar='is_domain', help='Whether to obtain TE domains, HiTE uses RepeatPeps.lib from RepeatMasker to obtain TE domains, 1: true, 0: false. default = [ ' + str(default_domain) + ' ]')
    parser.add_argument('--recover', metavar='is_recover', help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ ' + str(default_recover) + ' ]')
    parser.add_argument('--annotate', metavar='is_annotate', help='Whether to annotate the genome using the TE library generated, 1: true, 0: false. default = [ ' + str(default_annotate) + ' ]')
    parser.add_argument('--intact_anno', metavar='intact_anno', help='Whether to generate annotation of full-length TEs, 1: true, 0: false. default = [ ' + str(default_intact_annotate) + ' ]')
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
    parser.add_argument('--use_NeuralTE', metavar='use_NeuralTE', help='Whether to use NeuralTE to classify TEs, 1: true, 0: false. default = [' + str(default_use_NeuralTE) + ' ]')
    parser.add_argument('--is_wicker', metavar='is_wicker', help='Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [ ' + str(default_is_wicker) + ' ]')

    parser.add_argument('--flanking_len', metavar='flanking_len', help='The flanking length of candidates to find the true boundaries, default = [ ' + str(default_flanking_len) + ' ]')
    parser.add_argument('--fixed_extend_base_threshold', metavar='fixed_extend_base_threshold', help='The length of variation can be tolerated during pairwise alignment, default = [ '+str(default_fixed_extend_base_threshold)+' ]')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff', help='Cutoff of the candidates regarded as tandem region, default = [ '+str(default_tandem_region_cutoff)+' ]')
    parser.add_argument('--max_repeat_len', metavar='max_repeat_len', help='The maximum length of a single repeat, default = [ ' + str(default_max_single_repeat_len) + ' ]')
    parser.add_argument('--chrom_seg_length', metavar='chrom_seg_length', help='The length of genome segments, default = [ ' + str(default_chrom_seg_length) + ' ]')

    args = parser.parse_args()

    reference = args.genome
    threads = args.thread
    output_dir = args.outdir
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
    intact_anno = args.intact_anno
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
    use_NeuralTE = args.use_NeuralTE
    is_wicker = args.is_wicker

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
        output_dir = project_dir + '/output'
        log.logger.warning('\noutput directory path is empty, set to: ' + str(output_dir))

    if not os.path.isabs(reference):
        reference = os.path.abspath(reference)
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    if not os.path.exists(reference):
        log.logger.error('\nCannot find input genome assembly: ' + str(reference))
        parser.print_help()
        exit(-1)

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

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

    # if classified is None:
    #     classified = default_classified
    # else:
    #     classified = int(classified)

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

    if intact_anno is None:
        intact_anno = default_intact_annotate
    else:
        intact_anno = int(intact_anno)

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

    os.system('cp ' + reference + ' ' + tmp_output_dir)
    reference = tmp_output_dir + '/' + ref_filename

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

    LTR_harvest_parallel_Home = project_dir + '/bin/LTR_HARVEST_parallel'
    LTR_finder_parallel_Home = project_dir + '/bin/LTR_FINDER_parallel-master'
    NeuralTE_home = project_dir + '/bin/NeuralTE'
    EAHelitron = project_dir + '/bin/EAHelitron-master'
    HSDIR = project_dir + '/bin/HelitronScanner/TrainingSet'
    HSJAR = project_dir + '/bin/HelitronScanner/HelitronScanner.jar'
    rm2_script = project_dir + '/bin/get_family_summary_paper.sh'
    rm2_strict_script = project_dir + '/bin/get_family_summary_paper_0.99.sh'
    member_script_path = tools_dir + '/make_fasta_from_blast.sh'
    subset_script_path = tools_dir + '/ready_for_MSA.sh'
    sh_dir = project_dir + '/bin'
    lib_module = project_dir + '/library'
    protein_lib_path = project_dir + '/library/RepeatPeps.lib'

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
                    # '  [Setting] Is classified = [ ' + str(classified) + ' ] Default( ' + str(default_classified) + ' )\n'
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
                    '  [Setting] intact_anno = [ ' + str(intact_anno) + ' ]  Default( ' + str(default_intact_annotate) + ' )\n'                                                                                        
                    '  [Setting] search_struct = [ ' + str(search_struct) + ' ] Default( ' + str(default_search_struct) + ' )\n'                                                                                        
                    '  [Setting] BM_RM2 = [ ' + str(BM_RM2) + ' ]  Default( ' + str(default_BM_RM2) + ' )\n'
                    '  [Setting] BM_EDTA = [ ' + str(BM_EDTA) + ' ]  Default( ' + str(default_BM_EDTA) + ' )\n'
                    '  [Setting] BM_HiTE = [ ' + str(BM_HiTE) + ' ]  Default( ' + str(default_BM_HiTE) + ' )\n'
                    '  [Setting] EDTA_home = [' + str(EDTA_home) + ']\n'
                    '  [Setting] coverage_threshold = [ ' + str(coverage_threshold) + ' ]  Default( ' + str(default_coverage_threshold) + ' )\n'
                    '  [Setting] skip_HiTE = [ ' + str(skip_HiTE) + ' ]  Default( ' + str(default_skip_HiTE) + ' )\n'
                    '  [Setting] is_denovo_nonltr = [ ' + str(is_denovo_nonltr) + ' ]  Default( ' + str(default_is_denovo_nonltr) + ' )\n'
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
    library_dir = project_dir + '/library'

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

    pipeline_starttime = time.time()
    if skip_HiTE != 1:
        log.logger.info('Start step0: Structural Based LTR Searching')
        confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
        resut_file = confident_ltr_cut_path
        if not is_recover or not file_exist(resut_file):
            if te_type == 'all' or te_type == 'ltr':
                starttime = time.time()
                TEClass_home = project_dir + '/classification'
                LTR_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_LTR_transposons.py ' \
                                             + ' -g ' + reference + ' --ltrharvest_home ' + LTR_harvest_parallel_Home \
                                             + ' --ltrfinder_home ' + LTR_finder_parallel_Home + ' -t ' + str(threads) \
                                             + ' --tmp_output_dir ' + tmp_output_dir + ' --recover ' + str(recover) \
                                             + ' --miu ' + str(miu) + ' --use_NeuralTE ' + str(use_NeuralTE) + ' --is_wicker ' + str(is_wicker) \
                                             + ' --NeuralTE_home ' + NeuralTE_home + ' --TEClass_home ' + str(TEClass_home)
                log.logger.info(LTR_identification_command)
                os.system(LTR_identification_command)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step0: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        confident_other_path = tmp_output_dir + '/confident_other.fa'
        resut_file = confident_other_path
        if not is_recover or not file_exist(resut_file):
            if te_type == 'all' or te_type == 'non-ltr':
                starttime = time.time()
                log.logger.info('Start step1: homology-based other TE searching')
                other_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Other_transposons.py ' \
                                               + ' -r ' + reference \
                                               + ' -t ' + str(threads) \
                                               + ' --tmp_output_dir ' + tmp_output_dir  \
                                               + ' --library_dir ' + str(library_dir) + ' --recover ' + str(recover)
                log.logger.info(other_identification_command)
                os.system(other_identification_command)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step1: %.8s s" % (dtime))
        else:
            log.logger.info(resut_file + ' exists, skip...')

        # --------------------------------------------------------------------------------------
        starttime = time.time()
        log.logger.info('Start step2.0: Splitting genome assembly into chunks')
        split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
                                     + reference + ' --tmp_output_dir ' + tmp_output_dir \
                                     + ' --chrom_seg_length ' + str(chrom_seg_length) + ' --chunk_size ' + str(chunk_size)
        log.logger.info(split_genome_command)
        os.system(split_genome_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step2.0: %.8s s" % (dtime))

        reg_str = 'genome.cut(\d+).fa$'
        cut_references = []
        for filename in os.listdir(tmp_output_dir):
            match = re.search(reg_str, filename)
            if match:
                ref_index = match.group(1)
                cut_references.append((ref_index, tmp_output_dir + '/' + filename))

        # Using identified TEs to mask the genome in order to reduce computational load in all-vs-all alignments.
        prev_TE = tmp_output_dir + '/prev_TE.fa'
        if os.path.exists(confident_ltr_cut_path):
            os.system('cat ' + confident_ltr_cut_path + ' > ' + prev_TE)
        if curated_lib is not None and os.path.exists(curated_lib):
            os.system('cat ' + curated_lib + ' >> ' + prev_TE)
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
                    log.logger.info('Start 2.1: Coarse-grained boundary mapping')
                    coarse_boundary_command = 'cd ' + test_home + ' && python3 ' + test_home + '/coarse_boundary.py ' \
                                           + ' -g ' + cut_reference + ' --tmp_output_dir ' + tmp_output_dir \
                                           + ' --prev_TE ' + str(prev_TE) \
                                           + ' --fixed_extend_base_threshold ' + str(fixed_extend_base_threshold) \
                                           + ' --max_repeat_len ' + str(max_repeat_len) \
                                           + ' --thread ' + str(threads) \
                                           + ' --flanking_len ' + str(flanking_len) \
                                           + ' --tandem_region_cutoff ' + str(tandem_region_cutoff) \
                                           + ' --ref_index ' + str(ref_index) \
                                           + ' -r ' + reference + ' --recover ' + str(recover) \
                                           + ' --debug ' + str(debug)
                    log.logger.info(coarse_boundary_command)
                    os.system(coarse_boundary_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step2.1: %.8s s" % (dtime))
            else:
                log.logger.info(resut_file + ' exists, skip...')

            longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
            resut_file = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
            if not is_recover or not file_exist(resut_file):
                if te_type == 'all' or te_type == 'tir':
                    starttime = time.time()
                    log.logger.info('Start step2.2: determine fine-grained TIR')
                    tir_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_TIR_transposons.py ' \
                                                 + ' --seqs ' + longest_repeats_flanked_path \
                                                 + ' -t ' + str(threads)+' --TRsearch_dir ' + TRsearch_dir \
                                                 + ' --tmp_output_dir ' + tmp_output_dir \
                                                 + ' --tandem_region_cutoff ' + str(tandem_region_cutoff) \
                                                 + ' --ref_index ' + str(ref_index) \
                                                 + ' --subset_script_path ' + str(subset_script_path) \
                                                 + ' --plant ' + str(plant) \
                                                 + ' --flanking_len ' + str(flanking_len) \
                                                 + ' --recover ' + str(recover) \
                                                 + ' --debug ' + str(debug) \
                                                 + ' -r ' + reference \
                                                 + ' --split_ref_dir ' + split_ref_dir \
                                                 + ' --prev_TE  ' + prev_TE
                    log.logger.debug(tir_identification_command)
                    os.system(tir_identification_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step2.2: %.8s s" % (dtime))
            else:
                log.logger.info(resut_file + ' exists, skip...')

            resut_file = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
            if not is_recover or not file_exist(resut_file):
                if te_type == 'all' or te_type == 'helitron':
                    starttime = time.time()
                    log.logger.info('Start step2.3: determine fine-grained Helitron')
                    helitron_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Helitron_transposons.py --seqs ' \
                                                      + longest_repeats_flanked_path + ' -r ' + reference + ' -t ' + str(threads) \
                                                      + ' --tmp_output_dir ' + tmp_output_dir + ' --HSDIR ' + HSDIR + ' --HSJAR ' + HSJAR \
                                                      + ' --sh_dir ' + sh_dir + ' --EAHelitron ' + EAHelitron \
                                                      + ' --subset_script_path ' + subset_script_path \
                                                      + ' --ref_index ' + str(ref_index) + ' --flanking_len ' + str(flanking_len) \
                                                      + ' --recover ' + str(recover) + ' --debug ' + str(debug) + ' --split_ref_dir ' + split_ref_dir \
                                                      + ' --prev_TE  ' + prev_TE

                    log.logger.info(helitron_identification_command)
                    os.system(helitron_identification_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step2.3: %.8s s" % (dtime))
            else:
                log.logger.info(resut_file + ' exists, skip...')

            resut_file = tmp_output_dir + '/confident_non_ltr_' + str(ref_index) + '.fa'
            if not is_recover or not file_exist(resut_file):
                if te_type == 'all' or te_type == 'non-ltr':
                    starttime = time.time()
                    log.logger.info('Start step2.4: determine fine-grained Non-LTR')
                    non_ltr_identification_command = 'cd ' + test_home + ' && python3 ' + test_home + '/judge_Non_LTR_transposons.py'\
                                                 + ' --seqs ' + longest_repeats_flanked_path + ' -t ' + str(threads) \
                                                 + ' --subset_script_path ' + str(subset_script_path) \
                                                 + ' --tmp_output_dir ' + tmp_output_dir \
                                                 + ' --library_dir ' + str(library_dir) \
                                                 + ' --recover ' + str(recover) \
                                                 + ' --plant ' + str(plant) \
                                                 + ' --debug ' + str(debug) \
                                                 + ' --flanking_len ' + str(flanking_len) \
                                                 + ' --ref_index ' + str(ref_index) \
                                                 + ' --is_denovo_nonltr ' + str(is_denovo_nonltr) \
                                                 + ' -r ' + reference \
                                                 + ' --prev_TE  ' + prev_TE
                    log.logger.debug(non_ltr_identification_command)
                    os.system(non_ltr_identification_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step2.4: %.8s s" % (dtime))
            else:
                log.logger.info(resut_file + ' exists, skip...')

        # Merging all TEs generated from individual chunks.
        confident_tir_path = tmp_output_dir + '/confident_tir_merge.fa'
        confident_helitron_path = tmp_output_dir + '/confident_helitron_merge.fa'
        confident_non_ltr_path = tmp_output_dir + '/confident_non_ltr_merge.fa'
        os.system('rm -f ' + confident_tir_path)
        os.system('rm -f ' + confident_helitron_path)
        os.system('rm -f ' + confident_non_ltr_path)
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

        starttime = time.time()
        log.logger.info('Start step3: generate non-redundant library')
        TEClass_home = project_dir + '/classification'
        generate_lib_command = 'cd ' + test_home + ' && python3 ' + test_home + '/get_nonRedundant_lib.py' \
                               + ' --confident_ltr_cut ' + confident_ltr_cut_path \
                               + ' --confident_tir ' + confident_tir_path \
                               + ' --confident_helitron ' + confident_helitron_path \
                               + ' --confident_non_ltr ' + confident_non_ltr_path \
                               + ' --confident_other ' + confident_other_path \
                               + ' -t ' + str(threads) + ' --tmp_output_dir ' + tmp_output_dir \
                               + ' --test_home ' + str(test_home) + ' --use_NeuralTE ' + str(use_NeuralTE) \
                               + ' --is_wicker ' + str(is_wicker) \
                               + ' --NeuralTE_home ' + NeuralTE_home + ' --TEClass_home ' + str(TEClass_home) \
                               + ' --domain ' + str(domain) + ' --protein_path ' + str(protein_lib_path) \
                               + ' --curated_lib ' + str(curated_lib)
        log.logger.info(generate_lib_command)
        os.system(generate_lib_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step3: %.8s s" % (dtime))

    if intact_anno == 1:
        starttime = time.time()
        log.logger.info('Start step4: get full-length TE annotation')
        chr_name_map = tmp_output_dir + '/chr_name.map'
        ltr_list = tmp_output_dir + '/genome.rename.fa.pass.list'
        confident_tir_path = tmp_output_dir + '/confident_tir.fa'
        confident_helitron_path = tmp_output_dir + '/confident_helitron.fa'
        confident_non_ltr_path = tmp_output_dir + '/confident_non_ltr.fa'
        confident_other_path = tmp_output_dir + '/confident_other.fa'
        classified_TE_path = tmp_output_dir + '/TE_merge_tmp.fa.classified'
        full_length_anno_command = 'cd ' + test_home + ' && python3 ' + test_home + '/get_full_length_annotation.py' \
                                   + ' -t ' + str(threads) + ' --ltr_list ' + ltr_list \
                                   + ' --tir_lib ' + str(confident_tir_path) \
                                   + ' --helitron_lib ' + confident_helitron_path \
                                   + ' --nonltr_lib ' + confident_non_ltr_path \
                                   + ' --other_lib ' + confident_other_path \
                                   + ' --chr_name_map ' + chr_name_map \
                                   + ' -r ' + reference \
                                   + ' --module_home ' + test_home \
                                   + ' --tmp_output_dir ' + tmp_output_dir \
                                   + ' --TRsearch_dir ' + TRsearch_dir \
                                   + ' --search_struct ' + str(search_struct) \
                                   + ' --classified_TE_path ' + str(classified_TE_path)
        log.logger.info(full_length_anno_command)
        os.system(full_length_anno_command)

        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step4: %.8s s" % (dtime))

    confident_TE_consensus = tmp_output_dir + '/confident_TE.cons.fa'
    if not os.path.exists(confident_TE_consensus):
        log.logger.error('Error, Cannot find TE path: ' + confident_TE_consensus)
    else:
        starttime = time.time()
        log.logger.info('Start step5: annotate genome')
        TEClass_home = project_dir + '/classification'
        annotate_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/annotate_genome.py' \
                                  + ' -t ' + str(threads) + ' --classified_TE_consensus ' + confident_TE_consensus \
                                  + ' --annotate ' + str(annotate) \
                                  + ' -r ' + reference \
                                  + ' --tmp_output_dir ' + tmp_output_dir
        log.logger.info(annotate_genome_command)
        os.system(annotate_genome_command)

        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step5: %.8s s" % (dtime))

        starttime = time.time()
        log.logger.info('Start step6: Start conduct benchmarking of RepeatModeler2, EDTA, and HiTE')
        benchmarking_command = 'cd ' + test_home + ' && python3 ' + test_home + '/benchmarking.py' \
                            + ' --tmp_output_dir ' + tmp_output_dir \
                            + ' --BM_RM2 ' + str(BM_RM2) + ' --BM_EDTA ' + str(BM_EDTA) + ' --BM_HiTE ' + str(BM_HiTE) \
                            + ' --coverage_threshold ' + str(coverage_threshold) + ' -t ' + str(threads) + ' --lib_module ' + str(lib_module) \
                            + ' --TE_lib ' + str(confident_TE_consensus) + ' --rm2_script ' + str(rm2_script) \
                            + ' --rm2_strict_script ' + str(rm2_strict_script) + ' -r ' + reference
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
        log.logger.info("Running time of step6: %.8s s" % (dtime))

    clean_lib_command = 'cd ' + test_home + ' && python3 ' + test_home + '/clean_lib.py' \
                           + ' --tmp_output_dir ' + tmp_output_dir \
                           + ' --debug ' + str(debug)

    os.system(clean_lib_command)

    pipeline_endtime = time.time()
    dtime = pipeline_endtime - pipeline_starttime
    log.logger.info("Running time of the whole pipeline: %.8s s" % (dtime))
