import argparse
import os
import sys



current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)


from configs import config
from src.Util import read_fasta, store_fasta


def get_LTR_seq_from_bed(ref_contigs, bed_path, ltr_path):
    LTR_seqs = {}
    with open(bed_path, 'r') as f_r:
        for i, line in enumerate(f_r):
            if i == 0:
                continue
            else:
                line = line.replace('\n', '')
                parts = line.split('\t')
                chr_name = parts[0]
                LTR_start = parts[2]
                LTR_end = parts[5]
                if LTR_start != 'NA':
                    lLTR_start = int(parts[2])
                    lLTR_end = int(parts[3])
                    lLTR_seq = ref_contigs[chr_name][lLTR_start - 1: lLTR_end]
                    lLTR_name = chr_name + '_' + str(lLTR_start) + '-' + str(lLTR_end) + '-lLTR' + '#LTR'
                    LTR_seqs[lLTR_name] = lLTR_seq
                if LTR_end != 'NA':
                    rLTR_start = int(parts[4])
                    rLTR_end = int(parts[5])
                    rLTR_seq = ref_contigs[chr_name][rLTR_start - 1: rLTR_end]
                    rLTR_name = chr_name + '_' + str(rLTR_start) + '-' + str(rLTR_end) + '-rLTR' + '#LTR'
                    LTR_seqs[rLTR_name] = rLTR_seq
                if LTR_start != 'NA' and LTR_end != 'NA':
                    LTR_int_seq = ref_contigs[chr_name][lLTR_end: rLTR_start - 1]
                    LTR_int_name = chr_name + '_' + str(lLTR_end) + '-' + str(rLTR_start) + '-int' + '#LTR'
                    LTR_seqs[LTR_int_name] = LTR_int_seq
    f_r.close()
    store_fasta(LTR_seqs, ltr_path)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='generate Look4LTRs library')
    parser.add_argument('-i', metavar='Rtr_path', required=True,
                        help='The path of Rtr')
    parser.add_argument('-o', metavar='output_dir', required=True,
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('-g', required=True, metavar='genome', help='Input genome assembly path')



    args = parser.parse_args()
    genome_path = args.g
    output_dir = args.o
    Rtr_path = args.i

    project_dir = config.project_dir
    src_dir = project_dir + '/src'
    tool_dir = project_dir + '/tools'

    ref_names, ref_contigs = read_fasta(genome_path)

    # Step4. 将 scn 中的 left LTR提取出来
    ltr_path = output_dir + '/LTR.fa'
    ltr_cons = output_dir + '/LTR.cons'
    get_LTR_seq_from_bed(ref_contigs, Rtr_path, ltr_path)

    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_path + ' -o ' + ltr_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)