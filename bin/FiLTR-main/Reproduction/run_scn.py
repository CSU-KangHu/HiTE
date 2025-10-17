import argparse
import os
import sys
import time



current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from src.Util import convert_LtrDetector_scn, store_fasta, read_fasta
from configs import config


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

if __name__ == '__main__':
    tool_name = 'SCN TO library '
    version_num = '1.0.0'
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'

    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--genome', required=True, metavar='genome', help='Input genome assembly path')
    parser.add_argument('--scn_file', required=True, metavar='scn_file', help='Input scn_file')
    parser.add_argument('--threads', required=True, metavar='threads', help='threads')
    parser.add_argument('--out_dir', required=True, metavar='output_dir',
                        help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')

    args = parser.parse_args()
    genome_path = args.genome
    scn_file = args.scn_file
    threads = args.threads
    output_dir = args.out_dir

    ref_names, ref_contigs = read_fasta(genome_path)

    ltr_path = output_dir + '/LTR.fa'
    ltr_cons = output_dir + '/LTR.cons'
    get_LTR_seq_from_scn(ref_contigs, scn_file, ltr_path)

    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_path + ' -o ' + ltr_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)