import os
import time
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)


from configs import config
from src.Util import read_fasta, store_fasta, get_full_length_copies, getReverseSequence, generate_msa, read_fasta_v1

from concurrent.futures import ProcessPoolExecutor, as_completed

def get_both_ends_frame(candidate_sequence_path, ref_contigs, threads, temp_dir, output_dir, full_length_output_dir, split_ref_dir):
    debug = 0
    flanking_len = 100
    starttime = time.time()
    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(full_length_output_dir):
        os.makedirs(full_length_output_dir)

    # We are considering that the current running time is too long, maybe it is related to submitting one sequence for Blastn alignment at a time.
    # We will try to combine 10 sequences together and run Blastn once.
    # To increase CPU utilization, we will submit one thread to process 10 sequences.
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
    batch_member_files = []
    new_all_copies = {}
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
            copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
            if direct == '-':
                copy_seq = getReverseSequence(copy_seq)
            if len(copy_seq) < 100:
                continue
            new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end) + '(' + direct + ')'
            if not new_all_copies.__contains__(query_name):
                new_all_copies[query_name] = {}
            copy_contigs = new_all_copies[query_name]
            copy_contigs[new_name] = copy_seq
            new_all_copies[query_name] = copy_contigs
    for query_name in new_all_copies.keys():
        copy_contigs = new_all_copies[query_name]
        cur_member_file = temp_dir + '/' + query_name + '.blast.bed.fa'
        store_fasta(copy_contigs, cur_member_file)
        batch_member_files.append((query_name, cur_member_file))

    subset_script_path = os.getcwd() + '/tools/ready_for_MSA.sh'
    # Determine whether the multiple sequence alignment of each copied file satisfies the homology rule
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for batch_member_file in batch_member_files:
        job = ex.submit(generate_msa, batch_member_file, temp_dir, output_dir, full_length_output_dir, flanking_len, debug)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        left_frame_path, full_length_align_file = job.result()

    endtime = time.time()
    dtime = endtime - starttime

if __name__ == '__main__':
    # 提取repbase中所有的具有基因组的ltr元素
    species_genomes = {}
    with open('genome.info.bak', 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            species = parts[0]
            genome_path = parts[1]
            species_genomes[species] = genome_path

    all_repbase_path = '/public/home/hpc194701009/left_LTR_real_dataset/all_repbase.ref'
    repbase_names, repbase_contigs = read_fasta_v1(all_repbase_path)
    keep_ltr_contigs = {}
    total_ltr_count = 0
    ltr_tags = ['Gypsy', 'Copia', 'LTR Retrotransposon', 'BEL', 'LTR', 'Endogenous Retrovirus', 'Caulimoviridae']
    for name in repbase_names:
        parts = name.split('\t')
        if len(parts) < 3:
            continue
        raw_name = parts[0]
        label = parts[1]
        species = parts[2]
        if species in species_genomes and ('-LTR' in raw_name or '_LTR' in raw_name) and label in ltr_tags:
            total_ltr_count += 1
            if species not in keep_ltr_contigs:
                keep_ltr_contigs[species] = {}
            cur_species_ltr_contigs = keep_ltr_contigs[species]
            cur_species_ltr_contigs[raw_name] = repbase_contigs[name]
    print(total_ltr_count)

    repbase_ltr_dir = '/public/home/hpc194701009/left_LTR_real_dataset/repbase_ltr'
    if not os.path.exists(repbase_ltr_dir):
        os.makedirs(repbase_ltr_dir)
    for species in keep_ltr_contigs.keys():
        cur_species_ltr_contigs = keep_ltr_contigs[species]
        cur_repbase_ltr = repbase_ltr_dir + '/' + species + '.fa'
        store_fasta(cur_species_ltr_contigs, cur_repbase_ltr)


    # 本脚本的目的：
    # 1. 提取 Repbase 的LTR序列，并获取 LTR 的 左侧 100bp 多序列比对
    # 2. 提取 LtrDetector 的 FP 结果，并获取 LTR 的 左侧 100bp 多序列比对
    chrom_seg_length = 100000
    chunk_size = 40000
    threads = 40
    work_dir = '/public/home/hpc194701009/left_LTR_real_dataset'
    repbase_ltr_dir = work_dir + '/repbase_ltr'
    tmp_genome_dir = work_dir + '/genome'
    for species in species_genomes.keys():
        genome_path = species_genomes[species]
        genome_path = genome_path.replace('/public/home/hpc194701009', '/public/data/hpc174601028/hukang')
        cur_repbase_ltr = repbase_ltr_dir + '/' + species + '.fa'

        if not os.path.exists(tmp_genome_dir):
            os.makedirs(tmp_genome_dir)
        genome_name = os.path.basename(genome_path)
        os.system('cp ' + genome_path + ' ' + tmp_genome_dir)
        genome_path = tmp_genome_dir + '/' + genome_name

        if os.path.exists(genome_path) and os.path.exists(cur_repbase_ltr):
            tmp_output_dir = work_dir + '/positive/' + str(species).replace(' ', '_')
            if not os.path.exists(tmp_output_dir):
                os.makedirs(tmp_output_dir)

            genome_split_dir = tmp_output_dir + '/genome_chunks'
            if not os.path.exists(genome_split_dir):
                os.makedirs(genome_split_dir)

            # Step1. Splitting genome assembly into chunks
            test_home = os.getcwd()
            starttime = time.time()
            split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
                                   + genome_path + ' --tmp_output_dir ' + genome_split_dir \
                                   + ' --chrom_seg_length ' + str(chrom_seg_length) + ' --chunk_size ' + str(chunk_size)
            os.system(split_genome_command)
            endtime = time.time()
            dtime = endtime - starttime

            split_ref_dir = genome_split_dir + '/ref_chr'
            ref_contigs = {}
            for name in os.listdir(split_ref_dir):
                if name.endswith('.fa'):
                    cur_genome = split_ref_dir + '/' + name
                    cur_ref_names, cur_ref_contigs = read_fasta(cur_genome)
                    ref_contigs.update(cur_ref_contigs)

            # Step2. 获取 repbase LTR 终端序列 的 左侧框
            temp_dir = tmp_output_dir + '/repbase_ltr'
            output_dir = tmp_output_dir + '/positive'
            full_length_output_dir = tmp_output_dir + '/full_length'
            # get_left_frame(cur_repbase_ltr, ref_contigs, threads, temp_dir, output_dir, split_ref_dir)
            get_both_ends_frame(cur_repbase_ltr, ref_contigs, threads, temp_dir, output_dir, full_length_output_dir, split_ref_dir)

            if os.path.exists(temp_dir):
                os.system('rm -rf ' + temp_dir)
            if os.path.exists(genome_split_dir):
                os.system('rm -rf ' + genome_split_dir)

        if os.path.exists(genome_path):
            os.remove(genome_path)


    # chrom_seg_length = 100000
    # chunk_size = 40000
    # threads = 40
    # tmp_output_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/LtrDetector/dmel/output/left_LTR_sample'
    # genome_path = '/home/hukang/LTR_Benchmarking/LTR_libraries/LtrDetector/dmel/genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fa'
    #
    # # repbase_path = '/home/hukang/LtrHomo/library/maize.ltr.ref'
    # split_ref_dir = tmp_output_dir + '/ref_chr'
    # if not os.path.exists(tmp_output_dir):
    #     os.makedirs(tmp_output_dir)
    #
    # # Step1. Splitting genome assembly into chunks
    # test_home = os.getcwd()
    # starttime = time.time()
    # split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
    #                        + genome_path + ' --tmp_output_dir ' + tmp_output_dir \
    #                        + ' --chrom_seg_length ' + str(chrom_seg_length) + ' --chunk_size ' + str(chunk_size)
    # os.system(split_genome_command)
    # endtime = time.time()
    # dtime = endtime - starttime
    #
    # # Step2. 获取 repbase LTR 终端序列
    #
    # repbase_names, repbase_contigs = read_fasta(repbase_path)
    #
    # ltr_contigs = {}
    # for name in repbase_names:
    #     raw_name = name.split('#')[0]
    #     seq = repbase_contigs[name]
    #     if raw_name.endswith('-LTR') or raw_name.endswith('_LTR'):
    #         ltr_contigs[raw_name] = seq
    # repbase_ltr_seq = tmp_output_dir + '/repbase_ltr.fa'
    # store_fasta(ltr_contigs, repbase_ltr_seq)
    #
    # # Step3. 获取 repbase LTR 终端序列 的 左侧框
    # temp_dir = tmp_output_dir + '/repbase_ltr'
    # output_dir = tmp_output_dir + '/positive'
    # get_left_frame(repbase_ltr_seq, genome_path, threads, temp_dir, output_dir)
    #
    # # FP_path = '/home/hukang/LTR_Benchmarking/LTR_libraries/LtrDetector/dmel/output/FP.blastn.out'
    # # ltr_path = '/home/hukang/LTR_Benchmarking/LTR_libraries/LtrDetector/dmel/output/LTR.fa'
    # # # Step4. 获取 Ltrdetector FP 终端序列 的 左侧框
    # # FP_ltr_seq = tmp_output_dir + '/FP_ltr.fa'
    # # ltr_names, ltr_contigs = read_fasta(ltr_path)
    # # FP_ltr_contigs = {}
    # # ltr_terminal_names = set()
    # # with open(FP_path, 'r') as f_r:
    # #     for line in f_r:
    # #         parts = line.split('\t')
    # #         raw_name = parts[0]
    # #         chr_start = parts[2]
    # #         chr_end = parts[3]
    # #         if raw_name.endswith('-lLTR') and (chr_start in raw_name or chr_end in raw_name):
    # #             cur_name = raw_name+'#LTR'
    # #             seq = ltr_contigs[cur_name]
    # #             FP_ltr_contigs[raw_name] = seq
    # # store_fasta(FP_ltr_contigs, FP_ltr_seq)
    # #
    # # # Step5. 获取 Ltrdetector FP LTR 终端序列 的 左侧框
    # # temp_dir = tmp_output_dir + '/FP_ltr'
    # # output_dir = tmp_output_dir + '/negative'
    # # get_left_frame(FP_ltr_seq, genome_path, threads, temp_dir, output_dir, split_ref_dir)
    #
    #
