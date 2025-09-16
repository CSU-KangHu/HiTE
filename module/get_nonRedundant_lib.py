#!/usr/bin/env python
import argparse
import os
import shutil
import sys
import uuid

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, rename_fasta, remove_ltr_from_tir, get_domain_info, \
    ReassignInconsistentLabels, file_exist, create_or_clear_directory, copy_files, clean_old_tmp_files_by_dir, \
    ReassignUnknownLabels


def get_nonRedundant_lib(tmp_output_dir, confident_tir_path, confident_helitron_path, confident_non_ltr_path,
                         confident_other_path, confident_ltr_cut_path, threads, use_NeuralTE, is_wicker, curated_lib,
                         domain, log):
    test_home = cur_dir + '/module'
    NeuralTE_home = cur_dir + '/bin/NeuralTE'
    TEClass_home = cur_dir + '/classification'
    protein_path = cur_dir + '/library/RepeatPeps.lib'

    rename_fasta(confident_tir_path, confident_tir_path, 'TIR')
    rename_fasta(confident_helitron_path, confident_helitron_path, 'Helitron')
    rename_fasta(confident_non_ltr_path, confident_non_ltr_path, 'Denovo_Non-LTR')

    final_confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    final_confident_helitron_path = tmp_output_dir + '/confident_helitron.fa'
    final_confident_non_ltr_path = tmp_output_dir + '/confident_non_ltr.fa'
    final_confident_other_path = tmp_output_dir + '/confident_other.fa'

    # clustering and rename
    confident_tir_cons = confident_tir_path + '.cons'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_tir_path + ' -o ' + confident_tir_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    rename_fasta(confident_tir_cons, final_confident_tir_path, 'TIR')

    confident_helitron_cons = confident_helitron_path + '.cons'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_helitron_path + ' -o ' + confident_helitron_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    rename_fasta(confident_helitron_cons, final_confident_helitron_path, 'Helitron')

    confident_non_ltr_cons = confident_non_ltr_path + '.cons'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_non_ltr_path + ' -o ' + confident_non_ltr_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    rename_fasta(confident_non_ltr_cons, final_confident_non_ltr_path, 'Denovo_Non_LTR')

    # Merge all TE types (TIR+Helitron+Non_LTR+Other)
    confident_TE_path = tmp_output_dir + '/TE_merge_tmp.fa'
    if file_exist(final_confident_tir_path):
        os.system('cat ' + final_confident_tir_path + ' > ' + confident_TE_path)
    if file_exist(final_confident_helitron_path):
        os.system('cat ' + final_confident_helitron_path + ' >> ' + confident_TE_path)
    if file_exist(final_confident_non_ltr_path):
        os.system('cat ' + final_confident_non_ltr_path + ' >> ' + confident_TE_path)

    # Remove LTRs consist of other TE elements
    remove_ltr_from_tir(confident_TE_path, confident_ltr_cut_path, threads, tmp_output_dir)

    ref_rename_path = tmp_output_dir + '/genome.rename.fa'
    classified_TE_path = confident_TE_path + '.classified'
    # ClassifyTE except for LTR (LTRs have been classified)
    if use_NeuralTE:
        # classify LTR using NeuralTE
        NeuralTE_output_dir = tmp_output_dir + '/NeuralTE_all'
        if not os.path.exists(NeuralTE_output_dir):
            os.makedirs(NeuralTE_output_dir)
        NeuralTE_command = 'python ' + NeuralTE_home + '/src/Classifier.py --data ' + confident_TE_path \
                           + ' --genome ' + ref_rename_path + ' --use_TSD 1 --model_path ' \
                           + NeuralTE_home + '/models/NeuralTE-TSDs_model.h5 --out_dir ' \
                           + NeuralTE_output_dir + ' --thread ' + str(threads) + ' --is_wicker ' + str(is_wicker)
        log.logger.debug(NeuralTE_command)
        os.system(NeuralTE_command + ' > /dev/null 2>&1')
        NeuralTE_TE_path = NeuralTE_output_dir + '/classified_TE.fa'
        if os.path.exists(NeuralTE_TE_path):
            shutil.copy(NeuralTE_TE_path, classified_TE_path)
    else:
        # classify LTR using RepeatClassifier
        sample_name = 'test'
        TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                          + ' --consensus ' + confident_TE_path + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(threads) + ' -o ' + tmp_output_dir
        log.logger.debug(TEClass_command)
        os.system(TEClass_command)

    # merge classified LTRs and TEs
    confident_TE_path = tmp_output_dir + '/confident_TE.fa'
    if file_exist(classified_TE_path):
        shutil.copy(classified_TE_path, confident_TE_path)
    if file_exist(confident_ltr_cut_path):
        os.system('cat ' + confident_ltr_cut_path + ' >> ' + confident_TE_path)
    if file_exist(confident_other_path):
        os.system('cat ' + confident_other_path + ' >> ' + confident_TE_path)
    if curated_lib is not None:
        curated_lib = os.path.realpath(curated_lib)
        if file_exist(curated_lib):
            os.system('cat ' + curated_lib + ' >> ' + confident_TE_path)

    ReassignUnknownLabels(confident_TE_path)

    # # # Reassign Inconsistent Classification Labels
    # # ReassignInconsistentLabels(confident_TE_path)
    #
    # # classify all Unknown TEs with RepeatClassifier
    # unknown_tes = tmp_output_dir + '/unknown_TE.fa'
    # unknown_te_contigs = {}
    # te_names, te_contigs = read_fasta(confident_TE_path)
    # raw_name2full_name = {}
    # for name in te_names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     raw_name2full_name[raw_name] = name
    #     if 'Unknown' in name:
    #         unknown_te_contigs[raw_name] = te_contigs[name]
    # store_fasta(unknown_te_contigs, unknown_tes)
    #
    # sample_name = 'test'
    # TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + unknown_tes + ' --genome 1' \
    #                   + ' --thread_num ' + str(threads) + ' --split_num ' + str(threads) + ' -o ' + tmp_output_dir
    # log.logger.debug(TEClass_command)
    # os.system(TEClass_command)
    # classified_unknown_TE_path = unknown_tes + '.classified'
    # known_te_names, known_te_contigs = read_fasta(classified_unknown_TE_path)
    # for known_te_name in known_te_names:
    #     parts = known_te_name.split('#')
    #     if len(parts) == 2:
    #         raw_name = parts[0]
    #         unknown_name = raw_name2full_name[raw_name]
    #         unknown_seq = te_contigs[unknown_name]
    #         if unknown_name in te_contigs:
    #             del te_contigs[unknown_name]
    #         te_contigs[known_te_name] = unknown_seq
    # store_fasta(te_contigs, confident_TE_path)

    confident_TE_consensus = tmp_output_dir + '/confident_TE.cons.fa'
    if file_exist(confident_TE_path):
        cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                         + ' -G 0 -g 1 -A 80 -i ' + confident_TE_path + ' -o ' + confident_TE_consensus + ' -T 0 -M 0'
        os.system(cd_hit_command + ' > /dev/null 2>&1')

    # get domain for TEs
    if domain is not None and int(domain) == 1:
        output_table = confident_TE_consensus + '.domain'
        temp_dir = tmp_output_dir + '/domain'
        get_domain_info(confident_TE_consensus, protein_path, output_table, threads, temp_dir)

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE generating non-redundant library...')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--confident_ltr_cut', metavar='confident_ltr_cut',
                        help='The ltr path')
    parser.add_argument('--confident_tir', metavar='confident_tir',
                        help='The tir path')
    parser.add_argument('--confident_helitron', metavar='confident_helitron',
                        help='The helitron path')
    parser.add_argument('--confident_non_ltr', metavar='confident_non_ltr',
                        help='The denovo non-ltr path')
    parser.add_argument('--confident_other', metavar='confident_other',
                        help='The homology non-ltr path')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--use_NeuralTE', metavar='use_NeuralTE',
                        help='Whether to use NeuralTE to classify TEs, 1: true, 0: false.')
    parser.add_argument('--domain', metavar='domain',
                        help='Whether to obtain TE domains, HiTE uses RepeatPeps.lib from RepeatMasker to obtain TE domains, 1: true, 0: false.')
    parser.add_argument('--protein_path', metavar='protein_path',
                        help='The path of protein domain')
    parser.add_argument('--curated_lib', metavar='curated_lib',
                        help='The path of curated library')
    parser.add_argument('--is_wicker', metavar='is_wicker',
                        help='Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker.')
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    args = parser.parse_args()

    threads = int(args.t)
    confident_ltr_cut_path = args.confident_ltr_cut
    confident_tir_path = args.confident_tir
    confident_helitron_path = args.confident_helitron
    confident_non_ltr_path = args.confident_non_ltr
    confident_other_path = args.confident_other
    tmp_output_dir = args.tmp_output_dir
    use_NeuralTE = int(args.use_NeuralTE)
    domain = args.domain
    curated_lib = args.curated_lib
    is_wicker = args.is_wicker
    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)

    confident_ltr_cut_path = os.path.realpath(confident_ltr_cut_path)
    confident_tir_path = os.path.realpath(confident_tir_path)
    confident_helitron_path = os.path.realpath(confident_helitron_path)
    confident_non_ltr_path = os.path.realpath(confident_non_ltr_path)
    confident_other_path = os.path.realpath(confident_other_path)

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()

    tmp_output_dir = os.path.abspath(tmp_output_dir)

    log = Logger(tmp_output_dir+'/HiTE_lib.log', level='debug')

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'get_nonRedundant_lib_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)

        get_nonRedundant_lib(temp_dir, confident_tir_path, confident_helitron_path, confident_non_ltr_path,
                             confident_other_path, confident_ltr_cut_path, threads, use_NeuralTE, is_wicker, curated_lib,
                             domain, log)

        # 计算完之后将结果拷贝回输出目录
        copy_files(temp_dir, tmp_output_dir)

    except Exception as e:
        # 如果出现异常，打印错误信息并删除临时目录
        print(f"An error occurred: {e}")
        # if os.path.exists(temp_dir):
        #     shutil.rmtree(temp_dir)
        raise  # 重新抛出异常，以便上层代码可以处理

    else:
        # 如果没有异常，删除临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

