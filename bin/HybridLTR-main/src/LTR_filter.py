import argparse
import json
import os
import shutil
import subprocess
import sys
import time



current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from multiprocessing import cpu_count
from Util import read_fasta, store_fasta, Logger, read_scn, store_scn, get_LTR_seq_from_scn, get_recombination_ltr, \
    generate_both_ends_frame_from_seq, get_low_copy_LTR, judge_ltr_from_both_ends_frame, file_exist, \
    get_high_copy_LTR, alter_deep_learning_results, filter_tir, filter_sine, filter_helitron, \
    filter_ltr_by_flanking_cluster, filter_ltr_by_copy_num, filter_single_copy_ltr, remove_dirty_LTR, \
    filter_ltr_by_flank_seq_v2, deredundant_for_LTR_v5, get_all_potential_ltr_lines, \
    generate_both_ends_frame_from_seq_minimap2
from configs import config
from utils.data_util import sort_matrix_dir


def process_chunk(chunk, chunk_id, tmp_output_dir, flanking_len, threads, log, recover, debug,
                  is_flank_homo_cluster, is_flank_homo_cluster_both, is_handle_low_copy, is_use_homo_rule,
                  is_use_deep_model, is_filter_TIR, is_filter_Helitron, is_filter_SINE, reference, split_ref_dir,
                  max_copy_num, coverage_threshold, project_dir, src_dir, left_LTR_contigs, tool_dir):
    """
    处理单个块，生成中间结果文件。
    """
    # 创建当前块的临时目录
    chunk_dir = os.path.join(tmp_output_dir, f"chunk_{chunk_id}")
    os.makedirs(chunk_dir, exist_ok=True)

    # 生成当前块的 LTR 框架文件
    chunk_left_ltr_path = os.path.join(chunk_dir, "left_ltr.fasta")
    with open(chunk_left_ltr_path, "w") as f:
        for header, seq in chunk:
            f.write(f">{header}\n{seq}\n")

    # Step 3: Filter out false positive sequences based on the flanking regions of the terminal sequences.
    log.logger.info(f'Processing chunk {chunk_id}: Filter out false positive sequences based on the flanking regions of the terminal sequences.')
    confident_tir = os.path.join(tmp_output_dir, 'tir_' + str(chunk_id) + '.fa')
    confident_helitron = os.path.join(tmp_output_dir, 'helitron_' + str(chunk_id) + '.fa')
    confident_non_ltr = os.path.join(tmp_output_dir, 'non_ltr_' + str(chunk_id) + '.fa')
    confident_msa_file = os.path.join(tmp_output_dir, 'msa_flank_'+str(chunk_id)+'.txt')
    result_file = confident_msa_file

    if not recover or not file_exist(result_file):
        if file_exist(result_file):
            os.system('rm -f ' + result_file)

        # Perform a multiple sequence alignment of the regions flanking the LTR terminal sequence copies.
        temp_dir = os.path.join(chunk_dir, 'candidate_ltr')
        output_dir = os.path.join(chunk_dir, 'ltr_both_frames')
        full_length_output_dir = os.path.join(chunk_dir, 'full_length_frames')

        if not recover or not file_exist(output_dir) or not file_exist(full_length_output_dir):
            log.logger.debug('Generate LTR frames')
            # generate_both_ends_frame_from_seq(chunk_left_ltr_path, reference, flanking_len, threads, temp_dir, output_dir,
            #                                  full_length_output_dir, split_ref_dir, max_copy_num, coverage_threshold)

            generate_both_ends_frame_from_seq_minimap2(chunk_left_ltr_path, reference, flanking_len, threads, temp_dir,
                                                       output_dir, full_length_output_dir, max_copy_num)
        else:
            log.logger.info(output_dir + ' exists, skip...')
            log.logger.info(full_length_output_dir + ' exists, skip...')

        if is_flank_homo_cluster:
            # 去掉某一侧侧翼窗口高度同源的，通常是 truncated 的终端或重复区。
            temp_dir = os.path.join(chunk_dir, 'temp_sort')
            keep_output_dir = output_dir + '_sort'
            if not recover or not file_exist(keep_output_dir):
                log.logger.debug('Homologous clustering LTR frames')
                input_num, keep_num = sort_matrix_dir(output_dir, keep_output_dir, temp_dir, threads)
                log.logger.debug('Input LTR num: ' + str(input_num) + ', keep LTR num: ' + str(keep_num))
            else:
                log.logger.info(keep_output_dir + ' exists, skip...')

            if not debug:
                os.system('rm -rf ' + output_dir)
                os.system('rm -rf ' + temp_dir)

            output_dir = keep_output_dir

        if is_flank_homo_cluster_both:
            # 过滤掉具有高度同源的连接后的两侧侧翼区域，通常是LTR插入到其他TE中转座导致两侧区域完全一致，或者干脆就是假阳性。
            temp_dir = os.path.join(chunk_dir, 'temp_flank_cluster')
            keep_output_dir = output_dir + '_keep'
            if not recover or not file_exist(keep_output_dir):
                log.logger.debug('Homologous clustering flanking sequences of LTR')
                filter_ltr_by_flanking_cluster(output_dir, keep_output_dir, temp_dir, threads, log)
            else:
                log.logger.info(keep_output_dir + ' exists, skip...')

            if not debug:
                os.system('rm -rf ' + output_dir)
                os.system('rm -rf ' + temp_dir)

            output_dir = keep_output_dir

        output_path = os.path.join(chunk_dir, 'is_LTR.txt')
        if file_exist(output_path):
            os.remove(output_path)

        if is_handle_low_copy:
            lc_output_path = os.path.join(chunk_dir, 'is_LTR_homo.lc.txt')
            result_file = lc_output_path
            if not recover or not file_exist(result_file):
                log.logger.debug('Step 3.1: Use a rule-based method to filter out low-copy LTRs.')
                low_copy_output_dir = os.path.join(chunk_dir, 'low_copy_frames')
                get_low_copy_LTR(output_dir, low_copy_output_dir, threads, copy_num_threshold=5)

                type = 'Low copy'
                judge_ltr_from_both_ends_frame(low_copy_output_dir, lc_output_path, threads, type, flanking_len, log)

                if not debug:
                    os.system('rm -rf ' + low_copy_output_dir)
            else:
                log.logger.info(result_file + ' exists, skip...')
            os.system('cat ' + lc_output_path + ' >> ' + output_path)

        high_copy_output_dir = os.path.join(chunk_dir, 'high_copy_frames')
        result_file = high_copy_output_dir
        if not recover or not file_exist(result_file):
            log.logger.debug('Copy high-copy LTR frames for deep learning predicting')
            get_high_copy_LTR(output_dir, high_copy_output_dir, threads, copy_num_threshold=5)
        else:
            log.logger.info(result_file + ' exists, skip...')

        if not debug:
            os.system('rm -rf ' + output_dir)

        hc_output_path = os.path.join(chunk_dir, 'is_LTR_homo.hc.txt')
        if is_use_homo_rule:
            result_file = hc_output_path
            if not recover or not file_exist(result_file):
                type = 'High copy'
                judge_ltr_from_both_ends_frame(high_copy_output_dir, hc_output_path, threads, type, flanking_len, log)
            else:
                log.logger.info(result_file + ' exists, skip...')

        dl_output_path = os.path.join(chunk_dir, 'is_LTR_deep.txt')
        if is_use_deep_model:
            result_file = dl_output_path
            if not recover or not file_exist(result_file):
                file_names = os.listdir(high_copy_output_dir)
                if len(file_names) > 0:
                    model_path = os.path.join(project_dir, 'models/checkpoint_epoch_14.pth')
                    feature_output_dir = os.path.join(chunk_dir, 'feature_output_dir')
                    img_features = os.path.join(feature_output_dir, 'img_features.pt')
                    freq_features = os.path.join(feature_output_dir, 'freq_features.pt')
                    seq_names = os.path.join(feature_output_dir, 'seq_names.txt')
                    classify_command = f'python {src_dir}/Deep_Learning/hybridLTR_deep_main.py --matrix_dir {high_copy_output_dir} --threads {threads} --feature_output_dir {feature_output_dir} --model_path {model_path} --img_features {img_features} --freq_features {freq_features} --seq_names {seq_names} --output_dir {chunk_dir} --batch_size 256 --threshold 0.9 --device cpu'
                    log.logger.debug(classify_command)
                    os.system(classify_command)
            else:
                log.logger.info(result_file + ' exists, skip...')

        alter_dl_output_path = os.path.join(chunk_dir, 'is_LTR_deep.alter.txt')
        alter_deep_learning_results(dl_output_path, hc_output_path, alter_dl_output_path, high_copy_output_dir, log)
        os.system('cat ' + alter_dl_output_path + ' >> ' + output_path)

        if not debug:
            os.system('rm -rf ' + high_copy_output_dir)

        if is_filter_TIR:
            tir_output_path = os.path.join(chunk_dir, 'is_LTR_tir.txt')
            result_file = tir_output_path
            if not recover or not file_exist(result_file):
                filter_tir(output_path, tir_output_path, confident_tir, full_length_output_dir, threads, left_LTR_contigs, chunk_dir, tool_dir, flanking_len, log, debug)
            else:
                log.logger.info(result_file + ' exists, skip...')
            output_path = tir_output_path

        if is_filter_Helitron:
            helitron_output_path = os.path.join(chunk_dir, 'is_LTR_helitron.txt')
            result_file = helitron_output_path
            if not recover or not file_exist(result_file):
                filter_helitron(output_path, helitron_output_path, confident_helitron, full_length_output_dir, threads, left_LTR_contigs, chunk_dir, project_dir, flanking_len, log, debug)
            else:
                log.logger.info(result_file + ' exists, skip...')
            output_path = helitron_output_path

        if is_filter_SINE:
            sine_output_path = os.path.join(chunk_dir, 'is_LTR_sine.txt')
            result_file = sine_output_path
            if not recover or not file_exist(result_file):
                filter_sine(output_path, sine_output_path, confident_non_ltr, full_length_output_dir, threads, left_LTR_contigs, chunk_dir, flanking_len, log, debug)
            else:
                log.logger.info(result_file + ' exists, skip...')
            output_path = sine_output_path

        if not debug:
            os.system('rm -rf ' + full_length_output_dir)

        os.system('cat ' + output_path + ' > ' + confident_msa_file)
    else:
        log.logger.info(result_file + ' exists, skip...')

    return confident_msa_file, confident_tir, confident_helitron, confident_non_ltr


if __name__ == '__main__':
    tool_name = 'HybridLTR'
    version_num = '0.0.1'
    default_threads = int(cpu_count())
    default_miu = str(1.3e-8)
    default_flanking_len = 100

    default_is_filter_single = 1
    default_is_process_dirty_internal = 1
    default_is_use_flank_align = 1
    default_is_use_flank_MSA = 1
    default_is_flank_homo_cluster = 1
    default_is_flank_homo_cluster_both = 1
    default_is_filter_tandem = 1
    default_is_filter_TIR = 1
    default_is_filter_Helitron = 1
    default_is_filter_SINE = 1
    default_is_handle_low_copy = 1
    default_is_use_homo_rule = 1
    default_is_use_deep_model = 1
    default_is_use_structure = 0
    default_is_clean_internal = 1
    default_is_remove_nested = 0
    default_is_deredundant = 1
    default_is_output_lib = 1
    default_recover = 0

    default_skip_detect = 0
    default_debug = 0
    default_BM_RM2 = 0
    default_BM_EDTA = 0
    default_EDTA_home = ''
    default_BM_HiTE = 0
    default_coverage_threshold = 0.95


    # 1.parse args
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--scn', required=True, metavar='scn', help='Input scn file generated by LTR_finder and LTR_harvest')
    parser.add_argument('--genome', required=True, metavar='genome', help='Input genome assembly path')
    parser.add_argument('--out_dir', required=True, metavar='output_dir',
                        help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')

    parser.add_argument('--thread', metavar='thread_num',
                        help='Input thread num, default = [ ' + str(default_threads) + ' ]')
    parser.add_argument('--miu', metavar='miu',
                        help='The neutral mutation rate (per bp per ya), default = [ ' + str(default_miu) + ' ]')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='The flanking length of candidates to find the true boundaries, default = [ ' + str(default_flanking_len) + ' ]')
    parser.add_argument('--skip_detect', metavar='skip_detect',
                        help='Whether to skip_HiTE, 1: true, 0: false. default = [ ' + str(default_skip_detect) + ' ]')
    parser.add_argument('--debug', metavar='is_debug',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ ' + str(default_debug) + ' ]')

    parser.add_argument('--is_filter_single', metavar='is_filter_single',
                        help='Whether to filter LTRs with full-length copy number <= 1. default = [ ' + str(default_is_filter_single) + ' ]')
    parser.add_argument('--is_process_dirty_internal', metavar='is_process_dirty_internal',
                        help='Whether to filter recombined LTRs. default = [ ' + str(default_is_process_dirty_internal) + ' ]')
    parser.add_argument('--is_use_flank_align', metavar='is_use_flank_align',
                        help='Whether to filter out false positives that can be aligned to the LTR terminal flanking sequences. default = [ ' + str(default_is_use_flank_align) + ' ]')
    parser.add_argument('--is_use_flank_MSA', metavar='is_use_flank_MSA',
                        help='Whether to use the flanking region multiple sequence alignment strategy for filtering. default = [ ' + str(default_is_use_flank_MSA) + ' ]')
    parser.add_argument('--is_flank_homo_cluster', metavar='is_flank_homo_cluster',
                        help='Whether to perform homology clustering on flanking windows. default = [ ' + str(default_is_flank_homo_cluster) + ' ]')
    parser.add_argument('--is_flank_homo_cluster_both', metavar='is_flank_homo_cluster_both',
                        help='Whether to perform homology clustering on flanking windows. default = [ ' + str(default_is_flank_homo_cluster_both) + ' ]')

    parser.add_argument('--is_filter_tandem', metavar='is_filter_tandem',
                        help='Whether to filter LTR terminals composed of tandem repeats. default = [ ' + str(default_is_filter_tandem) + ' ]')
    parser.add_argument('--is_filter_TIR', metavar='is_filter_TIR',
                        help='Whether to filter LTR terminals composed of terminal inverted repeats (TIR). default = [ ' + str(default_is_filter_TIR) + ' ]')
    parser.add_argument('--is_filter_Helitron', metavar='is_filter_Helitron',
                        help='Whether to filter LTR terminals composed of Helitron. default = [ ' + str(default_is_filter_Helitron) + ' ]')
    parser.add_argument('--is_filter_SINE', metavar='is_filter_SINE',
                        help='Whether to filter LTR terminals composed of SINE. default = [ ' + str(default_is_filter_SINE) + ' ]')
    parser.add_argument('--is_handle_low_copy', metavar='is_handle_low_copy',
                        help='Whether to handle low-copy LTR separately. default = [ ' + str(default_is_handle_low_copy) + ' ]')
    parser.add_argument('--is_use_homo_rule', metavar='is_use_homo_rule',
                        help='Whether to use homology-rule to filter false positives. default = [ ' + str(default_is_use_homo_rule) + ' ]')
    parser.add_argument('--is_use_deep_model', metavar='is_use_deep_model',
                        help='Whether to use deep learning model. default = [ ' + str(default_is_use_deep_model) + ' ]')
    parser.add_argument('--is_use_structure', metavar='is_use_structure',
                        help='Whether to use deep learning model. default = [ ' + str(default_is_use_structure) + ' ]')
    parser.add_argument('--is_clean_internal', metavar='is_clean_internal',
                        help='Whether to clean LTR internal sequences. default = [ ' + str(default_is_clean_internal) + ' ]')
    parser.add_argument('--is_remove_nested', metavar='is_remove_nested',
                        help='Whether to unpack nested LTRs. default = [ ' + str(default_is_remove_nested) + ' ]')
    parser.add_argument('--is_deredundant', metavar='is_deredundant',
                        help='Whether to generate LTR cons. default = [ ' + str(default_is_deredundant) + ' ]')
    parser.add_argument('--is_output_lib', metavar='is_output_lib',
                        help='Whether to output LTR library. default = [ ' + str(default_is_output_lib) + ' ]')
    parser.add_argument('--recover', metavar='is_recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ ' + str(default_recover) + ' ]')

    parser.add_argument('--BM_RM2', metavar='BM_RM2',
                        help='Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false. default = [ ' + str(default_BM_RM2) + ' ]')
    parser.add_argument('--BM_EDTA', metavar='BM_EDTA',
                        help='Whether to conduct benchmarking of EDTA, 1: true, 0: false. default = [ ' + str(default_BM_EDTA) + ' ]')
    parser.add_argument('--BM_HiTE', metavar='BM_HiTE',
                        help='Whether to conduct benchmarking of HiTE, 1: true, 0: false. default = [ ' + str(default_BM_HiTE) + ' ]')
    parser.add_argument('--EDTA_home', metavar='EDTA_home',
                        help='When conducting benchmarking of EDTA, you will be asked to input EDTA home path.')
    parser.add_argument('--coverage_threshold', metavar='coverage_threshold',
                        help='The coverage threshold of benchmarking methods. default = [ ' + str(default_coverage_threshold) + ' ]')
    parser.add_argument('--species', metavar='species',
                        help='Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).')

    args = parser.parse_args()

    scn_file = args.scn
    reference = args.genome
    output_dir = args.out_dir
    threads = args.thread
    miu = args.miu
    flanking_len = args.flanking_len
    skip_detect = args.skip_detect
    debug = args.debug

    is_filter_single = args.is_filter_single
    is_process_dirty_internal = args.is_process_dirty_internal
    is_use_flank_align = args.is_use_flank_align
    is_use_flank_MSA = args.is_use_flank_MSA
    is_flank_homo_cluster = args.is_flank_homo_cluster
    is_flank_homo_cluster_both = args.is_flank_homo_cluster_both
    is_filter_tandem = args.is_filter_tandem
    is_filter_TIR = args.is_filter_TIR
    is_filter_Helitron = args.is_filter_Helitron
    is_filter_SINE = args.is_filter_SINE
    is_handle_low_copy = args.is_handle_low_copy
    is_use_homo_rule = args.is_use_homo_rule
    is_use_deep_model = args.is_use_deep_model
    is_use_structure = args.is_use_structure
    is_clean_internal = args.is_clean_internal
    is_remove_nested = args.is_remove_nested
    is_deredundant = args.is_deredundant
    is_output_lib = args.is_output_lib
    recover = args.recover

    BM_RM2 = args.BM_RM2
    BM_EDTA = args.BM_EDTA
    BM_HiTE = args.BM_HiTE
    EDTA_home = args.EDTA_home
    coverage_threshold = args.coverage_threshold
    species = args.species

    tmp_output_dir = os.path.abspath(output_dir + '/')
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    project_dir = config.project_dir
    src_dir = project_dir + '/src'
    tool_dir = project_dir + '/tools'

    log = Logger(tmp_output_dir + '/FiLTR.log', level='debug')

    if reference is None:
        log.logger.error('\nGenome path can not be empty')
        parser.print_help()
        exit(-1)
    if output_dir is None:
        output_dir = project_dir + '/output'
        log.logger.warning('\noutput directory path is empty, set to: ' + str(output_dir))

    if not os.path.isabs(scn_file):
        scn_file = os.path.abspath(scn_file)
    if not os.path.isabs(reference):
        reference = os.path.abspath(reference)
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    if flanking_len is None:
        flanking_len = default_flanking_len
    else:
        flanking_len = int(flanking_len)

    if miu is None:
        miu = default_miu
    else:
        miu = str(miu)

    if is_filter_single is None:
        is_filter_single = default_is_filter_single
    else:
        is_filter_single = int(is_filter_single)

    if is_process_dirty_internal is None:
        is_process_dirty_internal = default_is_process_dirty_internal
    else:
        is_process_dirty_internal = int(is_process_dirty_internal)

    if is_use_flank_align is None:
        is_use_flank_align = default_is_use_flank_align
    else:
        is_use_flank_align = int(is_use_flank_align)

    if is_use_flank_MSA is None:
        is_use_flank_MSA = default_is_use_flank_MSA
    else:
        is_use_flank_MSA = int(is_use_flank_MSA)

    if is_flank_homo_cluster is None:
        is_flank_homo_cluster = default_is_flank_homo_cluster
    else:
        is_flank_homo_cluster = int(is_flank_homo_cluster)

    if is_flank_homo_cluster_both is None:
        is_flank_homo_cluster_both = default_is_flank_homo_cluster_both
    else:
        is_flank_homo_cluster_both = int(is_flank_homo_cluster_both)

    if is_filter_tandem is None:
        is_filter_tandem = default_is_filter_tandem
    else:
        is_filter_tandem = int(is_filter_tandem)

    if is_filter_TIR is None:
        is_filter_TIR = default_is_filter_TIR
    else:
        is_filter_TIR = int(is_filter_TIR)

    if is_filter_Helitron is None:
        is_filter_Helitron = default_is_filter_Helitron
    else:
        is_filter_Helitron = int(is_filter_Helitron)

    if is_filter_SINE is None:
        is_filter_SINE = default_is_filter_SINE
    else:
        is_filter_SINE = int(is_filter_SINE)

    if is_handle_low_copy is None:
        is_handle_low_copy = default_is_handle_low_copy
    else:
        is_handle_low_copy = int(is_handle_low_copy)

    if is_use_homo_rule is None:
        is_use_homo_rule = default_is_use_homo_rule
    else:
        is_use_homo_rule = int(is_use_homo_rule)

    if is_use_deep_model is None:
        is_use_deep_model = default_is_use_deep_model
    else:
        is_use_deep_model = int(is_use_deep_model)

    if is_use_structure is None:
        is_use_structure = default_is_use_structure
    else:
        is_use_structure = int(is_use_structure)

    if is_clean_internal is None:
        is_clean_internal = default_is_clean_internal
    else:
        is_clean_internal = int(is_clean_internal)

    if is_remove_nested is None:
        is_remove_nested = default_is_remove_nested
    else:
        is_remove_nested = int(is_remove_nested)

    if is_deredundant is None:
        is_deredundant = default_is_deredundant
    else:
        is_deredundant = int(is_deredundant)

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

    if coverage_threshold is None:
        coverage_threshold = default_coverage_threshold
    else:
        coverage_threshold = float(coverage_threshold)

    if skip_detect is None:
        skip_HiTE = default_skip_detect
    else:
        skip_HiTE = int(skip_detect)

    if debug is None:
        debug = default_debug
    else:
        debug = int(debug)


    log.logger.info('Important Note: Please ensure that the chromosome names in the genome are in the format of Chr+number, such as Chr1 and Chr2.')

    split_ref_dir = tmp_output_dir + '/ref_chr'
    test_home = src_dir
    result_file = split_ref_dir
    if not recover or not file_exist(result_file):
        split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
                               + reference + ' --tmp_output_dir ' + tmp_output_dir
        log.logger.debug(split_genome_command)
        os.system(split_genome_command)
    else:
        log.logger.info(result_file + ' exists, skip...')

    ref_names, ref_contigs = read_fasta(reference)
    ltr_candidates, ltr_lines = read_scn(scn_file, log, remove_dup=True)

    dirty_dicts = {}
    # Step1: Filter records where the internal sequence of the current LTR contains its terminal sequences. This situation arises due to recombination involving shared terminals of two or more LTRs.
    if is_process_dirty_internal:
        remove_recomb_scn = tmp_output_dir + '/remove_recomb.scn'
        result_file = remove_recomb_scn
        if not recover or not file_exist(result_file):
            log.logger.info('Start step1: Process dirty internal LTR sequences')
            temp_dir = os.path.join(tmp_output_dir, 'get_recombination_ltr')
            os.makedirs(temp_dir, exist_ok=True)
            recombination_candidates = get_recombination_ltr(ltr_candidates, ref_contigs, threads, temp_dir, log)
            confident_lines = []
            for candidate_index in ltr_candidates.keys():
                if candidate_index not in recombination_candidates:
                    line = ltr_lines[candidate_index]
                    confident_lines.append(line)

            temp_path = tmp_output_dir + '/all_potential_ltr.json'
            temp_dir = os.path.join(tmp_output_dir, 'get_all_potential_ltr_lines')
            os.makedirs(temp_dir, exist_ok=True)
            confident_lines = get_all_potential_ltr_lines(confident_lines, reference, threads, temp_path, temp_dir, log)

            # 获取内部序列包含其他完整LTR的记录
            dirty_dicts = remove_dirty_LTR(confident_lines, log)

            store_scn(confident_lines, remove_recomb_scn)
            log.logger.debug('Current LTR num: ' + str(len(confident_lines)))

            if not debug:
                os.system('rm -f ' + temp_path)

        else:
            log.logger.info(result_file + ' exists, skip...')
        scn_file = remove_recomb_scn

    if is_use_flank_align:
        filter_terminal_align_scn = tmp_output_dir + '/filter_terminal_align.scn'
        result_file = filter_terminal_align_scn
        if not recover or not file_exist(result_file):
            log.logger.info('Start step2: Filter out false positives that can be aligned to the LTR terminal flanking sequences.')
            temp_dir = os.path.join(tmp_output_dir, 'filter_ltr_by_flank_seq_v2')
            os.makedirs(temp_dir, exist_ok=True)
            filter_ltr_by_flank_seq_v2(scn_file, filter_terminal_align_scn, reference, threads, temp_dir, log)
        else:
            log.logger.info(result_file + ' exists, skip...')
        scn_file = filter_terminal_align_scn

    # Based on the SCN file, extract the left LTR terminal sequences and create a mapping between the LTR sequences and the SCN lines.
    # This will facilitate the subsequent filtering of the SCN file.
    left_ltr_path = tmp_output_dir + '/left_LTR.fa'
    ltr_candidates, ltr_lines = read_scn(scn_file, log)
    leftLtr2Candidates = {}
    left_LTR_contigs = {}
    for candidate_index in ltr_candidates.keys():
        (chr_name, left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end) = ltr_candidates[candidate_index]
        if chr_name not in ref_contigs:
            log.logger.error('Error: Chromosome names in the SCN file do not match the input genome names. Please correct this and rerun.')
            exit(-1)
        ref_seq = ref_contigs[chr_name]
        left_ltr_name = chr_name + '-' + str(left_ltr_start) + '-' + str(left_ltr_end)
        left_ltr_seq = ref_seq[left_ltr_start-1: left_ltr_end]
        is_all_n = all(char == 'N' for char in left_ltr_seq)
        if is_all_n:
            continue
        # process duplicate name
        counter = 1
        original_name = left_ltr_name
        while left_ltr_name in left_LTR_contigs:
            left_ltr_name = f"{original_name}-rep{counter}"
            counter += 1
        left_LTR_contigs[left_ltr_name] = left_ltr_seq
        leftLtr2Candidates[left_ltr_name] = candidate_index
    store_fasta(left_LTR_contigs, left_ltr_path)

    if is_filter_tandem:
        left_ltr_filter_path = left_ltr_path + '.filter_tandem'
        result_file = left_ltr_filter_path
        if not recover or not file_exist(result_file):
            # Step 2: Filter out LTR termini primarily composed of tandem repeats.
            log.logger.info('Step 2: Filter out LTR termini primarily composed of tandem repeats.')
            tandem_filter_command = 'python ' + src_dir + '/filter_tandem_repeats.py -f ' + left_ltr_path + ' > ' + left_ltr_filter_path
            os.system(tandem_filter_command)
        else:
            log.logger.info(result_file + ' exists, skip...')

        if not debug:
            os.system('rm -f ' + left_ltr_path)

        left_ltr_path = left_ltr_filter_path
        filter_left_ltr_names, filter_left_ltr_contigs = read_fasta(left_ltr_path)
        log.logger.debug('Remove tandem LTR: ' + str(len(left_LTR_contigs) - len(filter_left_ltr_contigs)) + ', remaining LTR num: ' + str(len(filter_left_ltr_contigs)))

    # 将 left LTR 存成 scn 文件
    left_ltr_names, left_ltr_contigs = read_fasta(left_ltr_path)
    output_path = tmp_output_dir + '/left_ltr.scn'
    with open(output_path, 'w') as f_save:
        for cur_name in left_ltr_names:
            f_save.write(cur_name + '\t' + str(1) + '\n')

    if is_use_flank_MSA:
        # 将数据划分成多个块，防止中间文件过大导致磁盘空间不足
        max_copy_num = 100
        # 分块大小（每次处理的序列数）
        chunk_size = 5000
        result_files = []
        # 读取 left_ltr_path
        ltr_names, ltr_contigs = read_fasta(left_ltr_path)
        # 分块处理
        chunk = []
        chunk_id = 0
        for header in ltr_names:
            chunk.append((header, ltr_contigs[header]))
            if len(chunk) == chunk_size:
                chunk_result_file, confident_tir, confident_helitron, confident_non_ltr = process_chunk(chunk, chunk_id, tmp_output_dir, flanking_len, threads, log, recover,
                                                  debug, is_flank_homo_cluster, is_flank_homo_cluster_both,
                                                  is_handle_low_copy, is_use_homo_rule, is_use_deep_model,
                                                  is_filter_TIR, is_filter_Helitron, is_filter_SINE, reference,
                                                  split_ref_dir, max_copy_num, coverage_threshold, project_dir, src_dir,
                                                  left_LTR_contigs, tool_dir)
                result_files.append((chunk_result_file, confident_tir, confident_helitron, confident_non_ltr))

                # 清理当前块的中间文件
                if not debug:
                    chunk_dir = os.path.join(tmp_output_dir, f"chunk_{chunk_id}")
                    shutil.rmtree(chunk_dir)

                # 准备下一个块
                chunk = []
                chunk_id += 1

        # 处理最后一个块
        if chunk:
            chunk_result_file, confident_tir, confident_helitron, confident_non_ltr = process_chunk(chunk, chunk_id, tmp_output_dir, flanking_len, threads, log, recover,
                                              debug, is_flank_homo_cluster, is_flank_homo_cluster_both,
                                              is_handle_low_copy, is_use_homo_rule, is_use_deep_model, is_filter_TIR,
                                              is_filter_Helitron, is_filter_SINE, reference, split_ref_dir,
                                              max_copy_num, coverage_threshold, project_dir, src_dir, left_LTR_contigs,
                                              tool_dir)
            result_files.append((chunk_result_file, confident_tir, confident_helitron, confident_non_ltr))

            # 清理最后一个块的中间文件
            if not debug:
                chunk_dir = os.path.join(tmp_output_dir, f"chunk_{chunk_id}")
                shutil.rmtree(chunk_dir)

        # 合并所有块的结果文件
        output_path = os.path.join(tmp_output_dir, "msa_flank.txt")
        confident_tir_path = os.path.join(tmp_output_dir, "confident_tir_from_ltr.fa")
        confident_helitron_path = os.path.join(tmp_output_dir, "confident_helitron_from_ltr.fa")
        confident_non_ltr_path = os.path.join(tmp_output_dir, "confident_non_ltr_from_ltr.fa")
        total_result_files = [output_path, confident_tir_path, confident_helitron_path, confident_non_ltr_path]
        for i, cur_output_path in enumerate(total_result_files):
            with open(cur_output_path, "w") as outfile:
                for items in result_files:
                    item = items[i]
                    with open(item, "r") as infile:
                        shutil.copyfileobj(infile, outfile)

    if is_filter_single:
        intact_output_path = tmp_output_dir + '/intact_LTR_homo.txt'
        result_file = intact_output_path
        if not recover or not file_exist(intact_output_path):
            log.logger.info('Process LTRs with a full-length copy number <= 1')
            ltr_copies, internal_ltrs = filter_ltr_by_copy_num(output_path, leftLtr2Candidates, ltr_lines, reference, tmp_output_dir, split_ref_dir, threads, coverage_threshold, debug)

            if debug:
                # 把ltr_copies存成文件
                ltr_copies_json = tmp_output_dir + '/ltr_copies.json'
                with open(ltr_copies_json, 'w', encoding='utf-8') as f:
                    json.dump(ltr_copies, f, ensure_ascii=False, indent=4)

            # 处理单拷贝LTR。我们只保留具有 TSD 结构的单拷贝 且 内部具有完整protein的。
            lib_dir = project_dir + '/databases'
            ltr_protein_db = lib_dir + '/LTRPeps.lib'
            other_protein_db = lib_dir + '/OtherPeps.lib'
            single_copy_internals_file = tmp_output_dir + '/single_copy_internal.fa'
            filter_single_copy_ltr(intact_output_path, single_copy_internals_file, ltr_copies, internal_ltrs,
                                   ltr_protein_db, other_protein_db, tmp_output_dir, threads, reference,
                                   leftLtr2Candidates, ltr_lines, log, debug)
        else:
            log.logger.info(result_file + ' exists, skip...')
        output_path = intact_output_path

    # Step5. 过滤掉False positives
    FP_ltrs = {}
    true_ltrs = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            if not is_ltr:
                FP_ltrs[ltr_name] = is_ltr
            else:
                true_ltrs[ltr_name] = is_ltr

    confident_scn = tmp_output_dir + '/confident_ltr.scn'
    confident_lines = []
    for name in leftLtr2Candidates.keys():
        # if name not in FP_ltrs:
        if name in true_ltrs:
            candidate_index = leftLtr2Candidates[name]
            line = ltr_lines[candidate_index]
            confident_lines.append(line)
    store_scn(confident_lines, confident_scn)
    scn_file = confident_scn

    # Step6. 生成LTR library
    confident_ltr_terminal = tmp_output_dir + '/confident_ltr.terminal.fa'
    confident_ltr_internal = tmp_output_dir + '/confident_ltr.internal.fa'
    confident_ltr_intact = tmp_output_dir + '/intact_LTR.fa'
    ltr_intact_list = tmp_output_dir + '/intact_LTR.list'
    get_LTR_seq_from_scn(reference, scn_file, confident_ltr_terminal, confident_ltr_internal, confident_ltr_intact,
                         dirty_dicts, ltr_intact_list, miu)

    if is_output_lib:
        confident_ltr_path = tmp_output_dir + '/confident_ltr.fa'
        result_file = confident_ltr_path
        if not recover or not file_exist(result_file):
            confident_ltr_terminal_cons = confident_ltr_terminal + '.cons'
            terminal_coverage_threshold = 0.95
            if is_deredundant:
                # Step7. Remove redundancy from the LTR terminal results.
                starttime = time.time()
                log.logger.info('Start step6: Remove LTR terminal redundancy')
                type = 'terminal'
                # 对于 terminal 而言，其内部不太可能出现太大的插入删除变异，因此我们只是利用已识别的LTR终端序列去冗余，不再获取拷贝支持
                deredundant_for_LTR_v5(confident_ltr_terminal, tmp_output_dir, threads, type, terminal_coverage_threshold, debug)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step6: %.8s s" % (dtime))
            else:
                cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(terminal_coverage_threshold) \
                                 + ' -G 0 -g 1 -A 80 -i ' + confident_ltr_terminal + ' -o ' + confident_ltr_terminal_cons + ' -T 0 -M 0'
                os.system(cd_hit_command + ' > /dev/null 2>&1')
            # rename_fasta(confident_ltr_terminal_cons, confident_ltr_terminal_cons, 'LTR_terminal')

            if is_clean_internal:
                # Step6.2 清理内部序列。大量的内部序列是由于其他类型的TE插入，或者干脆内部序列就是大量串联重复，我们需要排除这种影响。
                log.logger.info('Start step7: Clean the internal sequences of LTRs by removing tandem repeats, LINEs, TIRs, and other elements (using proteins).')
                clean_internal_command = 'python ' + src_dir + '/clean_LTR_internal.py ' \
                                        + ' -t ' + str(threads) \
                                        + ' --tmp_output_dir ' + tmp_output_dir \
                                        + ' --internal_seq ' + confident_ltr_internal
                log.logger.debug(clean_internal_command)
                os.system(clean_internal_command)
                confident_ltr_internal += '.clean'

            confident_ltr_internal_cons = confident_ltr_internal + '.cons'
            internal_coverage_threshold = 0.8
            if is_deredundant:
                # Step7. Remove redundancy from the LTR internal results.
                starttime = time.time()
                log.logger.info('Start step6: Remove LTR internal redundancy')
                type = 'internal'
                # 对于 internal 而言，其内部可能出现较大的插入删除变异，因此我们要求比较宽松的0.8阈值
                deredundant_for_LTR_v5(confident_ltr_internal, tmp_output_dir, threads, type, internal_coverage_threshold, debug)
                endtime = time.time()
                dtime = endtime - starttime
                log.logger.info("Running time of step6: %.8s s" % (dtime))
            else:
                cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(internal_coverage_threshold) \
                                 + ' -G 0 -g 1 -A 80 -i ' + confident_ltr_internal + ' -o ' + confident_ltr_internal_cons + ' -T 0 -M 0'
                os.system(cd_hit_command + ' > /dev/null 2>&1')
            # rename_fasta(confident_ltr_internal_cons, confident_ltr_internal_cons, 'LTR_internal')

            # Step8. 生成一致性library
            os.system('cat ' + confident_ltr_terminal_cons + ' ' + confident_ltr_internal_cons + ' > ' + confident_ltr_path)

            if not debug:
                confident_ltr_terminal = tmp_output_dir + '/confident_ltr.terminal.fa'
                confident_ltr_internal = tmp_output_dir + '/confident_ltr.internal.fa'
                os.system('rm -f ' + confident_ltr_terminal + '.*')
                os.system('rm -f ' + confident_ltr_internal + '.*')
        else:
            log.logger.info(result_file + ' exists, skip...')

        # # Step9. 调用评估方法
        # evaluation_command = 'cd ' + project_dir + ' && python ' + src_dir + '/benchmarking.py --BM_RM2 ' + str(BM_RM2) + ' --BM_EDTA ' + str(BM_EDTA) + ' --BM_HiTE ' + str(BM_HiTE) + ' -t ' + \
        #                   str(threads) + ' --TE_lib ' + confident_ltr_path + ' -r ' + reference + \
        #                   ' --tmp_output_dir ' + tmp_output_dir + ' --recover ' + str(recover)
        # if EDTA_home is not None and EDTA_home != '':
        #     evaluation_command += ' --EDTA_home ' + EDTA_home
        # if species is not None:
        #     evaluation_command += ' --species ' + species
        # log.logger.debug(evaluation_command)
        # os.system(evaluation_command)

