import argparse
import os
from multiprocessing import cpu_count

from configs import config
from src.Util import Logger, read_fasta, split_chromosomes, store_fasta, file_exist, map_chr_position

current_folder = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.join(current_folder, ".")

if __name__ == '__main__':
    tool_name = 'HybridLTR_pipeline'
    version_num = '0.0.1'
    default_threads = int(cpu_count())
    default_miu = str(1.3e-8)
    default_is_output_lib = 1
    default_recover = 0

    default_skip_detect = 0
    default_debug = 0
    default_BM_RM2 = 0
    default_BM_EDTA = 0
    default_EDTA_home = ''
    default_BM_HiTE = 0

    # 1.parse args
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--genome', required=True, metavar='genome', help='Input genome assembly path')
    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')

    parser.add_argument('--thread', metavar='thread_num', help='Input thread num, default = [ ' + str(default_threads) + ' ]')
    parser.add_argument('--miu', metavar='miu', help='The neutral mutation rate (per bp per ya), default = [ ' + str(default_miu) + ' ]')
    parser.add_argument('--skip_detect', metavar='skip_detect', help='Whether to skip_HybridLTR, 1: true, 0: false. default = [ ' + str(default_skip_detect) + ' ]')
    parser.add_argument('--debug', metavar='is_debug', help='Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ ' + str(default_debug) + ' ]')
    parser.add_argument('--is_output_lib', metavar='is_output_lib', help='Whether to output LTR library. default = [ ' + str(default_is_output_lib) + ' ]')
    parser.add_argument('--recover', metavar='is_recover', help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ ' + str(default_recover) + ' ]')

    parser.add_argument('--BM_RM2', metavar='BM_RM2', help='Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false. default = [ ' + str(default_BM_RM2) + ' ]')
    parser.add_argument('--BM_EDTA', metavar='BM_EDTA', help='Whether to conduct benchmarking of EDTA, 1: true, 0: false. default = [ ' + str(default_BM_EDTA) + ' ]')
    parser.add_argument('--BM_HiTE', metavar='BM_HiTE', help='Whether to conduct benchmarking of HiTE, 1: true, 0: false. default = [ ' + str(default_BM_HiTE) + ' ]')
    parser.add_argument('--EDTA_home', metavar='EDTA_home', help='When conducting benchmarking of EDTA, you will be asked to input EDTA home path.')
    parser.add_argument('--species', metavar='species', help='Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).')

    args = parser.parse_args()

    reference = args.genome
    output_dir = args.out_dir
    threads = args.thread
    miu = args.miu
    skip_detect = args.skip_detect
    debug = args.debug
    is_output_lib = args.is_output_lib
    recover = args.recover

    BM_RM2 = args.BM_RM2
    BM_EDTA = args.BM_EDTA
    BM_HiTE = args.BM_HiTE
    EDTA_home = args.EDTA_home
    species = args.species

    tmp_output_dir = os.path.abspath(output_dir + '/')
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    project_dir = config.project_dir
    src_dir = project_dir + '/src'
    tool_dir = project_dir + '/tools'

    log = Logger(tmp_output_dir + '/HybridLTR_pipeline.log', level='debug')

    if reference is None:
        log.logger.error('\nGenome path can not be empty')
        parser.print_help()
        exit(-1)
    if output_dir is None:
        output_dir = project_dir + '/output'
        log.logger.warning('\noutput directory path is empty, set to: ' + str(output_dir))

    if not os.path.isabs(reference):
        reference = os.path.abspath(reference)
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    if miu is None:
        miu = default_miu
    else:
        miu = str(miu)

    if is_output_lib is None:
        is_output_lib = default_is_output_lib
    else:
        is_output_lib = int(is_output_lib)

    if recover is None:
        recover = default_recover
    else:
        recover = int(recover)

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

    if skip_detect is None:
        skip_HiTE = default_skip_detect
    else:
        skip_HiTE = int(skip_detect)

    if debug is None:
        debug = default_debug
    else:
        debug = int(debug)

    # Step 0. 将 reference 进行切割，确保最长的序列不超过 10 Mb
    chunk_size = 10_000_000
    chromosomes_names, chromosomes_dict = read_fasta(reference)
    chromosomes_dict, position_map = split_chromosomes(chromosomes_dict, max_length=chunk_size)
    split_reference = output_dir + '/genome_split.fa'
    store_fasta(chromosomes_dict, split_reference)

    scn_file = output_dir + '/total.scn'
    result_file = scn_file
    if not recover or not file_exist(result_file):
        # Step 1. 调用 LTR_detector 并行版本
        log.logger.debug('Start running LTR_detector to obtain candidate LTRs')
        LTR_detector_output_dir = output_dir + '/LTR_detector_output'
        LTR_detector_command = 'cd ' + project_dir + '/Reproduction && python run_LtrDetector_parallel.py --genome ' \
                               + split_reference + ' --genome_dir '+ LTR_detector_output_dir + ' --out_dir ' + LTR_detector_output_dir \
                               + ' --LtrDetector_home ' + project_dir + '/bin/LtrDetector/bin --threads ' + str(threads)
        log.logger.debug(LTR_detector_command)
        os.system(LTR_detector_command + ' > /dev/null 2>&1')

        LTR_detector_scn_file = LTR_detector_output_dir + '/total.scn'
        # 将分割的染色体名称映射回原始的名称和位置
        map_chr_position(LTR_detector_scn_file, scn_file, position_map)
        # if not debug:
        #     os.system('rm -rf ' + LTR_detector_output_dir)
    else:
        log.logger.info(result_file + ' exists, skip...')


    log.logger.debug('Start running HybridLTR to get confident LTRs')
    # Step 2. 调用 HybridLTR 进行 LTR 过滤
    HybridLTR_command = 'cd ' + project_dir + '/src && python LTR_filter.py --scn ' + scn_file + ' --genome ' \
                        + reference + ' --out_dir ' + output_dir + ' --thread ' + str(threads) \
                        + ' --recover ' + str(recover) + ' --debug ' + str(debug) + ' --BM_RM2 ' \
                        + str(BM_RM2) + ' --BM_HiTE ' + str(BM_HiTE)  + ' --BM_EDTA ' + str(BM_EDTA) \
                        + ' --debug ' + str(debug) + ' --miu ' + str(miu) + ' --is_output_lib ' + str(is_output_lib)
    if EDTA_home is not None and EDTA_home.strip() != '':
        HybridLTR_command += ' --EDTA_home ' + str(EDTA_home)
    if species is None or species.strip() == '':
        HybridLTR_command += ' --species test'
    else:
        HybridLTR_command += ' --species ' + str(species)
    log.logger.debug(HybridLTR_command)
    os.system(HybridLTR_command)
