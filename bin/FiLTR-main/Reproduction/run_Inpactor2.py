import argparse
import os
import sys
import time

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from src.Util import read_fasta, store_fasta, read_fasta_v1
from configs import config

def split_fasta(cur_path, output_dir, num_chunks):
    split_files = []

    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    os.makedirs(output_dir)

    names, contigs = read_fasta_v1(cur_path)
    num_names = len(names)
    chunk_size = num_names // num_chunks

    for i in range(num_chunks):
        chunk_start = i * chunk_size
        chunk_end = chunk_start + chunk_size if i < num_chunks - 1 else num_names
        chunk = names[chunk_start:chunk_end]
        output_path = output_dir + '/out_' + str(i) + '.fa'
        with open(output_path, 'w') as out_file:
            for name in chunk:
                seq = contigs[name]
                out_file.write('>'+name+'\n'+seq+'\n')
        split_files.append(output_path)
    return split_files

def run_command(command):
    subprocess.run(command, check=True, shell=True)

def parse_ltr_log(ltr_log):
    ltr_pos = {}
    cur_seq_name = None
    with open(ltr_log, 'r') as f_r:
        for line in f_r:
            if line.startswith('load sequence'):
                seq_name = line.split(' ')[-1].split('\t')[0]
            elif 'Length ltr=' in line:
                parts = line.strip().split(' ')[0].split('..')
                lLTR_info = parts[0].replace('(', '').replace(')', '').split(',')
                rLTR_info = parts[1].replace('(', '').replace(')', '').split(',')
                lLTR_start = lLTR_info[0]
                lLTR_end = lLTR_info[1]
                rLTR_start = rLTR_info[0]
                rLTR_end = rLTR_info[1]
                pos = (lLTR_start, lLTR_end, rLTR_start, rLTR_end)
                if seq_name is not None:
                    ltr_pos[seq_name] = pos
                    seq_name = None
    return ltr_pos

def identify_terminals(split_file, output_dir, tool_dir):
    ltr_file = split_file + '.ltr'
    ltr_log = ltr_file + '.log'
    tir_file = split_file + '.itr'
    tir_log = tir_file + '.log'

    ltrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/ltrsearch -l 100 ' + split_file + ' > ' + ltr_log
    itrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + split_file+ ' > ' + tir_log
    run_command(ltrsearch_command)
    run_command(itrsearch_command)
    # os.system(ltrsearch_command)
    # os.system(itrsearch_command)

    # 解析.log 文件，读取 'load sequence' 和 'Length ltr=' 标志
    LTR_info = parse_ltr_log(ltr_log)

    # Update the header of the split_file, adding two columns LTR:1-206,4552-4757 TIR:1-33,3869-3836.
    update_split_file = split_file + '.updated'
    update_contigs = {}
    names, contigs = read_fasta_v1(split_file)
    for name in names:
        orig_name = name.split('\t')[0]
        LTR_str = 'LTR:'
        lLTR_start, lLTR_end, rLTR_start, rLTR_end = LTR_info[orig_name]
        LTR_str += str(lLTR_start) + '-' + str(lLTR_end) + ',' + str(rLTR_start) + '-' + str(rLTR_end)
        update_name = name + '\t' + LTR_str
        update_contigs[update_name] = contigs[name]
    store_fasta(update_contigs, update_split_file)
    return update_split_file

def generate_terminal_info(data_path, work_dir, tool_dir, threads):
    output_dir = work_dir + '/temp'
    # Split the file into threads blocks.
    split_files = split_fasta(data_path, output_dir, threads)

    # Parallelize the identification of LTR and TIR.
    cur_update_path = data_path + '.update'
    os.system('rm -f ' + cur_update_path)
    with ProcessPoolExecutor(threads) as executor:
        futures = []
        for split_file in split_files:
            future = executor.submit(identify_terminals, split_file, output_dir, tool_dir)
            futures.append(future)
        executor.shutdown(wait=True)

        is_exit = False
        for future in as_completed(futures):
            update_split_file = future.result()
            if isinstance(update_split_file, str):
                os.system('cat ' + update_split_file + ' >> ' + cur_update_path)
            else:
                print(f"An error occurred: {update_split_file}")
                is_exit = True
                break
        if is_exit:
            print('Error occur, exit...')
            exit(1)
        else:
            os.system('mv ' + cur_update_path + ' ' + data_path)

    return data_path


if __name__ == '__main__':
    tool_name = 'Generating Inpactor2 library'
    version_num = '1.0.0'
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'

    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--inpactor2_result', required=True, metavar='inpactor2_result', help='Input inpactor2 result')
    parser.add_argument('--threads', required=True, metavar='threads', help='threads')
    parser.add_argument('--out_dir', required=True, metavar='output_dir',
                        help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')
    args = parser.parse_args()
    inpactor2_result = args.inpactor2_result
    threads = int(args.threads)
    output_dir = args.out_dir

    project_dir = config.project_dir
    src_dir = project_dir + '/src'
    tool_dir = project_dir + '/tools'

    work_dir = output_dir
    Inpactor2_fasta = inpactor2_result
    names, contigs = read_fasta(Inpactor2_fasta)
    # 去掉标签存储
    no_label_Inpactor2_fasta = Inpactor2_fasta + '.updated'
    new_contigs = {}
    for name in names:
        raw_name = name.split('#')[0]
        new_contigs[raw_name] = contigs[name]
    store_fasta(new_contigs, no_label_Inpactor2_fasta)
    data = generate_terminal_info(no_label_Inpactor2_fasta, work_dir, tool_dir, threads)

    # 根据 header 中的 LTR:1-979,6445-5482 分离出 terminal 和 internal
    ltr_path = work_dir + '/LTR.fa'
    ltr_cons = work_dir + '/LTR.cons'
    names, contigs = read_fasta_v1(data)
    LTR_seqs = {}
    for name in names:
        parts = name.split('\t')
        raw_name = parts[0]
        pos_parts = parts[1].split(':')[1].split(',')
        left_terminal_parts = pos_parts[0].split('-')
        right_terminal_parts = pos_parts[1].split('-')
        lLTR_start = int(left_terminal_parts[0])
        lLTR_end = int(left_terminal_parts[1])
        rLTR_start = int(right_terminal_parts[0])
        rLTR_end = int(right_terminal_parts[1])
        lLTR_seq = contigs[name][lLTR_start - 1: lLTR_end]
        rLTR_seq = contigs[name][rLTR_start - 1: rLTR_end]
        LTR_int_seq = contigs[name][lLTR_end: rLTR_start - 1]
        lLTR_name = raw_name + '-lLTR' + '#LTR'
        LTR_seqs[lLTR_name] = lLTR_seq
        LTR_int_name = raw_name + '-int' + '#LTR'
        LTR_seqs[LTR_int_name] = LTR_int_seq
        rLTR_name = raw_name + '-rLTR' + '#LTR'
        LTR_seqs[rLTR_name] = rLTR_seq
    store_fasta(LTR_seqs, ltr_path)
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_path + ' -o ' + ltr_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)