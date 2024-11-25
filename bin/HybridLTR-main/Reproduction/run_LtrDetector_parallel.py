import argparse
import os
import sys
import time



current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from src.Util import convert_LtrDetector_scn, store_fasta, read_fasta, split_chromosomes
from configs import config
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_LTR_seq_from_scn(ref_contigs, scn_path, ltr_path):
    LTR_seqs = {}
    with open(scn_path, 'r') as f_r:
        for i, line in enumerate(f_r):
            if line.startswith('#'):
                continue
            else:
                line = line.replace('\n', '')
                parts = line.split(' ')
                LTR_start = int(parts[0])
                LTR_end = int(parts[1])
                chr_name = parts[11]
                lLTR_start = int(parts[3])
                lLTR_end = int(parts[4])
                rLTR_start = int(parts[6])
                rLTR_end = int(parts[7])
                lLTR_seq = ref_contigs[chr_name][lLTR_start-1: lLTR_end]
                rLTR_seq = ref_contigs[chr_name][rLTR_start - 1: rLTR_end]
                LTR_int_seq = ref_contigs[chr_name][lLTR_end: rLTR_start - 1]

                # LTR_seq = ref_contigs[chr_name][LTR_start-1: LTR_end]
                # LTR_name = chr_name + '_' + str(LTR_start) + '-' + str(LTR_end) + '#LTR'
                lLTR_name = chr_name + '_' + str(LTR_start) + '-' + str(LTR_end) + '-lLTR' + '#LTR'
                LTR_seqs[lLTR_name] = lLTR_seq
                LTR_int_name = chr_name + '_' + str(LTR_start) + '-' + str(LTR_end) + '-int' + '#LTR'
                LTR_seqs[LTR_int_name] = LTR_int_seq
                rLTR_name = chr_name + '_' + str(LTR_start) + '-' + str(LTR_end) + '-rLTR' + '#LTR'
                LTR_seqs[rLTR_name] = rLTR_seq
    f_r.close()
    store_fasta(LTR_seqs, ltr_path)

def run_LtrDetector(fasta_dir, LtrDetector_home, output_dir):
    # Step2. run LtrDetector parallel
    LtrDetector_command = LtrDetector_home + '/LtrDetector -fasta ' + fasta_dir + ' -destDir ' + output_dir + ' -nThreads ' + str(1)
    os.system(LtrDetector_command)
    return 1

if __name__ == '__main__':
    tool_name = 'LtrDetector parallel'
    version_num = '1.0.0'
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'

    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--genome', required=True, metavar='genome', help='Input genome assembly path')
    parser.add_argument('--genome_dir', required=True, metavar='genome_dir', help='Input genome_dir')
    parser.add_argument('--threads', required=True, metavar='threads', help='threads')
    parser.add_argument('--out_dir', required=True, metavar='output_dir',
                        help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')
    parser.add_argument('--LtrDetector_home', required=True, metavar='LtrDetector_home',
                        help='The path of LtrDetector_home')

    args = parser.parse_args()
    genome_path = args.genome
    genome_dir = args.genome_dir
    threads = int(args.threads)
    output_dir = args.out_dir
    LtrDetector_home = args.LtrDetector_home

    chrom_seg_length = 100000
    chunk_size = 40000

    project_dir = config.project_dir
    src_dir = project_dir + '/src'
    tool_dir = project_dir + '/tools'

    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    os.makedirs(output_dir)


    # Step1. Splitting genome assembly into chunks
    test_home = current_folder
    starttime = time.time()
    split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
                           + genome_path + ' --tmp_output_dir ' + genome_dir
    os.system(split_genome_command)
    endtime = time.time()
    dtime = endtime - starttime

    split_ref_dir = genome_dir + '/ref_chr'

    # ref_contigs = {}
    # for name in os.listdir(split_ref_dir):
    #     if name.endswith('.fa'):
    #         cur_genome = split_ref_dir + '/' + name
    #         cur_ref_names, cur_ref_contigs = read_fasta(cur_genome)
    #         ref_contigs.update(cur_ref_contigs)

    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index in os.listdir(split_ref_dir):
        fasta_dir = split_ref_dir + '/' + str(partition_index)
        obj = ex.submit(run_LtrDetector, fasta_dir, LtrDetector_home, output_dir)
        objs.append(obj)
    ex.shutdown(wait=True)
    for obj in as_completed(objs):
        cur_res = obj.result()


    # Step3. 合并所有的LtrDetector 结果, 转换 LTRDetector 的输出为scn格式
    total_bed = output_dir + '/total.bed'
    os.system('cd ' + output_dir +  ' && cat *.bed > ' + total_bed)

    scn_file = output_dir + '/total.scn'
    convert_LtrDetector_scn(total_bed, scn_file)

    # # Step4. 将 scn 中的 left LTR提取出来
    # ltr_path = output_dir + '/LTR.fa'
    # ltr_cons = output_dir + '/LTR.cons'
    # get_LTR_seq_from_scn(ref_contigs, scn_file, ltr_path)
    #
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + ltr_path + ' -o ' + ltr_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)