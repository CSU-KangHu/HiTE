#!/usr/bin/env python
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import uuid
from datetime import datetime
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from concurrent.futures import ProcessPoolExecutor, as_completed
from Util import Logger, copy_files, create_or_clear_directory, \
    store_fasta, read_fasta, get_full_length_copies, getReverseSequence, run_find_members_v8, rename_fasta, \
    ReassignInconsistentLabels, file_exist, lib_add_prefix


def filter_detected_TEs(temp_dir, threads, low_copy_file, panTE_lib, TE_type, log):
    tmp_file = os.path.join(temp_dir, TE_type + '.merge.fa')
    tmp_cons = tmp_file + '.cons'
    os.system('cat ' + low_copy_file + ' ' + panTE_lib + ' > ' + tmp_file)
    cd_hit_command = 'cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -d 0 -g 1 -A 80 -i ' + tmp_file + ' -o ' + tmp_cons + ' -T ' + str(threads) + ' -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')

    keep_low_copy_file = os.path.join(temp_dir, TE_type + '.keep.fa')
    cons_names, cons_contigs = read_fasta(tmp_cons)
    # 解析聚类文件中的单拷贝
    cluster_file = tmp_cons + '.clstr'
    cluster_idx = -1
    clusters = {}
    with open(cluster_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('>'):
                cluster_idx = line.split(' ')[1]
            else:
                if not clusters.__contains__(cluster_idx):
                    clusters[cluster_idx] = []
                cur_cluster = clusters[cluster_idx]
                name = line.split(',')[1].split(' ')[1].strip()[1:]
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
    # 判断簇中是否包含 panTE 序列，例如header中包含#
    # 从一致性序列中去掉包含 panTE 簇中所包含的序列，剩下的序列即为需要恢复的低拷贝序列
    removed_ids = set()
    for cluster_idx in clusters.keys():
        cur_cluster = clusters[cluster_idx]
        is_panTE_clstr = False
        for seq_name in cur_cluster:
            if '#' in seq_name:
                is_panTE_clstr = True
                break
        if is_panTE_clstr:
            for seq_name in cur_cluster:
                removed_ids.add(seq_name)
    for seq_name in removed_ids:
        if seq_name in cons_contigs:
            del cons_contigs[seq_name]
    store_fasta(cons_contigs, keep_low_copy_file)
    return keep_low_copy_file

def for_test(output_dir, threads, panTE_lib, reference, genome_name, log):
    full_length_copies_path = f'{output_dir}/{genome_name}.full_length.copies'
    full_length_gff_path = f'{output_dir}/{genome_name}.full_length.gff'
    gff_path = f'{output_dir}/{genome_name}.gff'
    os.system('touch ' + full_length_copies_path)
    os.system('touch ' + full_length_gff_path)
    os.system('touch ' + gff_path)


def get_single_genome_copies(candidate_sequence_path, split_ref_dir, flanking_len, debug):
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(candidate_sequence_path)
    split_files = []
    cur_contigs = {}
    for i, name in enumerate(names):
        cur_file = temp_dir + '/' + str(batch_id) + '.fa'
        cur_contigs[name] = contigs[name]
        if len(cur_contigs) == batch_size:
            store_fasta(cur_contigs, cur_file)
            split_files.append(cur_file)
            cur_contigs = {}
            batch_id += 1
    if len(cur_contigs) > 0:
        cur_file = temp_dir + '/' + str(batch_id) + '.fa'
        store_fasta(cur_contigs, cur_file)
        split_files.append(cur_file)
        batch_id += 1

    ref_contigs = {}
    # 遍历目录下的所有文件
    for filename in os.listdir(split_ref_dir):
        if filename.endswith('.fa'):
            file_path = os.path.join(split_ref_dir, filename)
            cur_names, cur_contigs = read_fasta(file_path)
            ref_contigs.update(cur_contigs)

    # ref_names, ref_contigs = read_fasta(reference)
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies, cur_split_files, split_ref_dir, debug)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_all_copies = job.result()
        all_copies.update(cur_all_copies)
    # extend copies
    extend_all_copies = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        for copy in copies:
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            direct = copy[4]
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            extend_copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
            if direct == '-':
                extend_copy_seq = getReverseSequence(extend_copy_seq)
            if len(extend_copy_seq) < 100:
                continue
            new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end) + '(' + direct + ')'
            if not extend_all_copies.__contains__(query_name):
                extend_all_copies[query_name] = {}
            copy_contigs = extend_all_copies[query_name]
            copy_contigs[new_name] = extend_copy_seq
            extend_all_copies[query_name] = copy_contigs
    return extend_all_copies

def get_pan_genome_copies(keep_tir_low_copy, keep_helitron_low_copy, keep_non_ltr_low_copy, genome_paths, temp_dir):
    keep_tir_extend_copies = {}
    keep_helitron_extend_copies = {}
    keep_non_ltr_extend_copies = {}
    raw_tir_names, raw_tir_contigs = read_fasta(keep_tir_low_copy)
    raw_helitron_names, raw_helitron_contigs = read_fasta(keep_helitron_low_copy)
    raw_non_ltr_names, raw_non_ltr_contigs = read_fasta(keep_non_ltr_low_copy)
    for i, (genome_name, genome_path) in enumerate(genome_paths):
        cur_temp_dir = os.path.join(temp_dir, genome_name)
        create_or_clear_directory(cur_temp_dir)

        (ref_dir, ref_filename) = os.path.split(genome_path)
        reference = cur_temp_dir + '/' + ref_filename
        os.system('cp ' + genome_path + ' ' + reference)

        split_genome_command = 'split_genome_chunks.py -g ' \
                               + reference + ' --tmp_output_dir ' + cur_temp_dir \
                               + ' --chrom_seg_length ' + str(100000) + ' --chunk_size ' + str(400)
        log.logger.info(split_genome_command)
        os.system(split_genome_command)
        split_ref_dir = cur_temp_dir + '/ref_chr'

        flanking_len = 50
        debug = 0
        # 获取TIR拷贝
        tir_extend_copies = get_single_genome_copies(keep_tir_low_copy, split_ref_dir, flanking_len, debug)
        remain_tir_contigs = {}
        for query_name in tir_extend_copies.keys():
            copy_contigs = tir_extend_copies[query_name]
            # 如果当前获得了超过100个拷贝或者当前已经是最后一轮，则不需要进行下一轮
            if len(copy_contigs) >= 100 or i == len(genome_paths) - 1:
                keep_tir_extend_copies[query_name] = copy_contigs
            else:
                if query_name in raw_tir_contigs:
                    remain_tir_contigs[query_name] = raw_tir_contigs[query_name]
        # 经过一轮获取拷贝之后，剩下的序列仍然没有获得100个拷贝，因此需要进一步比对到其他基因组上
        store_fasta(remain_tir_contigs, keep_tir_low_copy)

        # 获取Helitron拷贝
        helitron_extend_copies = get_single_genome_copies(keep_helitron_low_copy, split_ref_dir, flanking_len, debug)
        remain_helitron_contigs = {}
        for query_name in helitron_extend_copies.keys():
            copy_contigs = helitron_extend_copies[query_name]
            # 如果当前获得了超过100个拷贝或者当前已经是最后一轮，则不需要进行下一轮
            if len(copy_contigs) >= 100 or i == len(genome_paths) - 1:
                keep_helitron_extend_copies[query_name] = copy_contigs
            else:
                if query_name in raw_helitron_contigs:
                    remain_helitron_contigs[query_name] = raw_helitron_contigs[query_name]
        # 经过一轮获取拷贝之后，剩下的序列仍然没有获得100个拷贝，因此需要进一步比对到其他基因组上
        store_fasta(remain_helitron_contigs, keep_helitron_low_copy)

        # 获取non_ltr拷贝
        non_ltr_extend_copies = get_single_genome_copies(keep_non_ltr_low_copy, split_ref_dir, flanking_len, debug)
        remain_non_ltr_contigs = {}
        for query_name in non_ltr_extend_copies.keys():
            copy_contigs = non_ltr_extend_copies[query_name]
            # 如果当前获得了超过100个拷贝或者当前已经是最后一轮，则不需要进行下一轮
            if len(copy_contigs) >= 100 or i == len(genome_paths) - 1:
                keep_non_ltr_extend_copies[query_name] = copy_contigs
            else:
                if query_name in raw_non_ltr_contigs:
                    remain_non_ltr_contigs[query_name] = raw_non_ltr_contigs[query_name]
        # 经过一轮获取拷贝之后，剩下的序列仍然没有获得100个拷贝，因此需要进一步比对到其他基因组上
        store_fasta(remain_non_ltr_contigs, keep_non_ltr_low_copy)

        # 记得删除目录，防止由于泛基因组数量过多/大，导致/tmp磁盘空间不足
        if os.path.exists(cur_temp_dir):
            shutil.rmtree(cur_temp_dir)

    # 将那些经过泛基因组之后，获得>=5拷贝存成文件
    tir_batch_member_files = []
    cur_tir_temp_dir = os.path.join(temp_dir, 'tir_members')
    os.makedirs(cur_tir_temp_dir, exist_ok=True)
    for query_name in keep_tir_extend_copies.keys():
        copy_contigs = keep_tir_extend_copies[query_name]
        if len(copy_contigs) >= 5:
            valid_query_filename = re.sub(r'[<>:"/\\|?*]', '-', query_name)
            extend_member_file = cur_tir_temp_dir + '/' + valid_query_filename + '.blast.bed.fa'
            store_fasta(copy_contigs, extend_member_file)
            query_seq = raw_tir_contigs[query_name]
            tir_batch_member_files.append((query_name, query_seq, extend_member_file))

    helitron_batch_member_files = []
    cur_helitron_temp_dir = os.path.join(temp_dir, 'helitron_members')
    os.makedirs(cur_helitron_temp_dir, exist_ok=True)
    for query_name in keep_helitron_extend_copies.keys():
        copy_contigs = keep_helitron_extend_copies[query_name]
        if len(copy_contigs) >= 5:
            valid_query_filename = re.sub(r'[<>:"/\\|?*]', '-', query_name)
            extend_member_file = cur_helitron_temp_dir + '/' + valid_query_filename + '.blast.bed.fa'
            store_fasta(copy_contigs, extend_member_file)
            query_seq = raw_helitron_contigs[query_name]
            helitron_batch_member_files.append((query_name, query_seq, extend_member_file))

    non_ltr_batch_member_files = []
    cur_non_ltr_temp_dir = os.path.join(temp_dir, 'non_ltr_members')
    os.makedirs(cur_non_ltr_temp_dir, exist_ok=True)
    for query_name in keep_non_ltr_extend_copies.keys():
        copy_contigs = keep_non_ltr_extend_copies[query_name]
        if len(copy_contigs) >= 5:
            valid_query_filename = re.sub(r'[<>:"/\\|?*]', '-', query_name)
            extend_member_file = cur_non_ltr_temp_dir + '/' + valid_query_filename + '.blast.bed.fa'
            store_fasta(copy_contigs, extend_member_file)
            query_seq = raw_non_ltr_contigs[query_name]
            non_ltr_batch_member_files.append((query_name, query_seq, extend_member_file))

    pan_genome_copies_dict = {}
    pan_genome_copies_dict['tir'] = (cur_tir_temp_dir, tir_batch_member_files)
    pan_genome_copies_dict['helitron'] = (cur_helitron_temp_dir, helitron_batch_member_files)
    pan_genome_copies_dict['non_ltr'] = (cur_non_ltr_temp_dir, non_ltr_batch_member_files)
    return pan_genome_copies_dict


def filter_true_TEs(batch_member_files, real_TEs, temp_dir, subset_script_path):
    # Determine whether the multiple sequence alignment of each copied file satisfies the homology rule
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for batch_member_file in batch_member_files:
        job = ex.submit(run_find_members_v8, batch_member_file, temp_dir, subset_script_path,
                        1, TE_type, 0, result_type='cons')
        jobs.append(job)
    ex.shutdown(wait=True)

    not_found_boundary = 0
    full_length1 = 0
    copy_nums = {}
    true_te_names = set()
    true_tes = {}
    for job in as_completed(jobs):
        result_info = job.result()
        cur_name, cur_seq, info, copy_count, extend_member_file = result_info
        if info == 'nb':
            not_found_boundary += 1
        elif info == 'fl1':
            full_length1 += 1
        elif info.startswith('copy_num:'):
            copy_num = int(info.split('copy_num:')[1])
            if not copy_nums.__contains__(copy_num):
                copy_nums[copy_num] = 0
            cur_copy_num = copy_nums[copy_num]
            cur_copy_num += 1
            copy_nums[copy_num] = cur_copy_num
        if cur_name is not None:
            true_tes[cur_name] = cur_seq
            true_te_names.add(cur_name)
    store_fasta(true_tes, real_TEs)


if __name__ == "__main__":
    default_use_NeuralTE = 1

    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE recover low copy TEs.")
    parser.add_argument("--genome_name", type=str, help="Name of the genome.")
    parser.add_argument("--tir_low_copy", type=str, help="low copy tir path, to recover tir using pan-genome")
    parser.add_argument("--helitron_low_copy", type=str, help="low copy helitron path, to recover helitron using pan-genome")
    parser.add_argument("--non_ltr_low_copy", type=str, help="low copy non_ltr path, to recover non_ltr using pan-genome")
    parser.add_argument("--panTE_lib", type=str, help="pan TE lib.")
    parser.add_argument("--genome_list", type=str, help="Path to the genome list file.")
    parser.add_argument("--pan_genomes_dir", type=str, help="Path to the directory containing pan-genomes.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")
    parser.add_argument('--use_NeuralTE', type=int, default=1,
                        help='Whether to use NeuralTE to classify TEs, 1: true, 0: false.')
    parser.add_argument('--is_wicker', type=int, default=0,
                        help='Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker.')

    # 解析参数
    args = parser.parse_args()
    raw_genome_name = args.genome_name
    tir_low_copy = args.tir_low_copy
    helitron_low_copy = args.helitron_low_copy
    non_ltr_low_copy = args.non_ltr_low_copy
    panTE_lib = args.panTE_lib
    genome_list = args.genome_list
    pan_genomes_dir = args.pan_genomes_dir
    threads = args.threads
    use_NeuralTE = args.use_NeuralTE
    is_wicker = args.is_wicker

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/pan_recover_low_copy_TEs.log', level='debug')

    # 1.1. 读取 genomes
    genome_paths = []
    with open(genome_list, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            genome_name = parts[0]
            cur_genome_path = os.path.join(pan_genomes_dir, genome_name)
            if not os.path.isabs(cur_genome_path):
                cur_genome_path = os.path.abspath(cur_genome_path)
            if not os.path.exists(cur_genome_path):
                log.logger.error(f'Cannot find genome path: {cur_genome_path}')
                sys.exit(-1)
            genome_paths.append((genome_name, cur_genome_path))


    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = '/tmp/pan_recover_low_copy_TEs_' + str(unique_id)
    create_or_clear_directory(temp_dir)
    # 利用 cd-hit-est 去掉低拷贝中包含在 panTE library 中的序列
    TE_type = 'tir'
    keep_tir_low_copy = filter_detected_TEs(temp_dir, threads, tir_low_copy, panTE_lib, TE_type, log)
    TE_type = 'helitron'
    keep_helitron_low_copy = filter_detected_TEs(temp_dir, threads, helitron_low_copy, panTE_lib, TE_type, log)
    TE_type = 'non_ltr'
    keep_non_ltr_low_copy = filter_detected_TEs(temp_dir, threads, non_ltr_low_copy, panTE_lib, TE_type, log)
    # 将保留下来的低拷贝候选TE放入一个队列，将队列中的序列依次比对回泛基因组获取拷贝，当低拷贝候选TE 的拷贝数超过 100 个时，将该低拷贝TE移出队列，继续剩余的TE比对，直至队列为空。
    # 函数返回一个dict: key 为 type, value为 batch_member_files 列表（记录每条序列对应的全长拷贝）
    get_copies_dir = os.path.join(temp_dir, 'get_copies')
    pan_genome_copies_dict = get_pan_genome_copies(keep_tir_low_copy, keep_helitron_low_copy, keep_non_ltr_low_copy, genome_paths, get_copies_dir)

    subset_script_path = project_dir + '/tools/ready_for_MSA.sh'
    # 依次对tir, helitron, non_ltr 类型的TE调用动态边界调整算法，过滤出真实的TE
    recover_panTE_lib = os.path.join(temp_dir, 'panTE.recover.fa')
    for TE_type in pan_genome_copies_dict.keys():
        cur_temp_dir, batch_member_files = pan_genome_copies_dict[TE_type]
        real_TEs = os.path.join(temp_dir, 'real_' + TE_type + '.fa')
        filter_true_TEs(batch_member_files, real_TEs, cur_temp_dir, subset_script_path)

        confident_recover_TE_path = os.path.join(temp_dir, 'confident_recover_' + TE_type + '.fa')
        # 生成非冗余序列
        real_TEs_cons = real_TEs + '.cons'
        cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                         + ' -G 0 -g 1 -A 80 -i ' + real_TEs + ' -o ' + real_TEs_cons + ' -T 0 -M 0' + ' > /dev/null 2>&1'
        os.system(cd_hit_command)

        prefix_map = {'tir': 'TIR_Recover', 'helitron': 'Helitron_Recover', 'non_ltr': 'Non-LTR_Recover'}
        rename_fasta(real_TEs_cons, confident_recover_TE_path, prefix_map[TE_type])
        os.system('cat ' + confident_recover_TE_path + ' >> ' + recover_panTE_lib)

    # 调用分类器对恢复的TE进行分类
    NeuralTE_home = project_dir + '/bin/NeuralTE'
    TEClass_home = project_dir + '/classification'
    recover_classified_lib = recover_panTE_lib + '.classified'
    if use_NeuralTE:
        # classify LTR using NeuralTE
        NeuralTE_output_dir = temp_dir + '/NeuralTE_recover'
        os.makedirs(NeuralTE_output_dir, exist_ok=True)
        NeuralTE_command = 'python ' + NeuralTE_home + '/src/Classifier.py --data ' + recover_panTE_lib \
                           + ' --model_path ' + NeuralTE_home + '/models/NeuralTE_model.h5 --out_dir ' \
                           + NeuralTE_output_dir + ' --thread ' + str(threads) + ' --is_wicker ' + str(is_wicker)
        log.logger.debug(NeuralTE_command)
        os.system(NeuralTE_command + ' > /dev/null 2>&1')
        NeuralTE_TE_path = NeuralTE_output_dir + '/classified_TE.fa'
        os.system('cp ' + NeuralTE_TE_path + ' ' + recover_classified_lib)
    else:
        # classify LTR using RepeatClassifier
        sample_name = 'test'
        TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                          + ' --consensus ' + recover_panTE_lib + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(threads) + ' -o ' + temp_dir
        log.logger.debug(TEClass_command)
        os.system(TEClass_command)

    # Reassign Inconsistent Classification Labels
    ReassignInconsistentLabels(recover_classified_lib)

    raw_name = raw_genome_name.split('.')[0]
    # 为文件加前缀
    if file_exist(recover_classified_lib):
        lib_add_prefix(recover_classified_lib, raw_name)

    # 计算完之后将结果拷贝回输出目录
    # copy_files(temp_dir, output_dir)
    shutil.copytree(temp_dir, output_dir, dirs_exist_ok=True)

    # 删除临时目录
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)