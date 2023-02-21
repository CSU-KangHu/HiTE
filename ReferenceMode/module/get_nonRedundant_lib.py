import argparse
import os
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, rename_fasta


def remove_self_alignment(blastnResults_path, new_blastnResults_path):
    records = []
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            else:
                records.append(line)
    with open(new_blastnResults_path, 'w') as f_save:
        for line in records:
            f_save.write(line)
    f_save.close()

def remove_nest(blastnResults_path, query_path, output, coverage = 0.95, identity_threshold = 95):

    minlen = 100
    new_query_contigs = {}
    query_names, query_contigs = read_fasta(query_path)
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            line = str(line).replace('\n', '')
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            if not query_contigs.__contains__(query_name) or not query_contigs.__contains__(subject_name):
                continue
            query_len = len(query_contigs[query_name])
            subject_len = len(query_contigs[subject_name])
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if s_start > s_end:
                tmp = s_start
                s_start = s_end
                s_end = tmp
            if query_name == subject_name:
                continue
            if float(alignment_len)/query_len < coverage or identity < identity_threshold:
                continue
            subject_seq = query_contigs[subject_name]
            sbj_seq_p1 = subject_seq[0: s_start]
            sbj_seq_p2 = subject_seq[s_end:]
            sbj_seq_new = sbj_seq_p1 + sbj_seq_p2
            sbj_len_new = len(sbj_seq_new)
            if sbj_len_new >= minlen:
                if not new_query_contigs.__contains__(subject_name):
                    new_query_contigs[subject_name] = sbj_seq_new
                else:
                    old_sbj = new_query_contigs[subject_name]
                    if sbj_len_new < len(old_sbj):
                        new_query_contigs[subject_name] = sbj_seq_new

    for query_name in query_contigs.keys():
        if new_query_contigs.__contains__(query_name):
            query_contigs[query_name] = new_query_contigs[query_name]
    store_fasta(query_contigs, output)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--sample_name', metavar='sample_name',
                        help='e.g., rice')
    parser.add_argument('--blast_program_dir', metavar='blast_program_dir',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2')
    parser.add_argument('--RepeatModeler_Home', metavar='RepeatModeler_Home',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/RepeatModeler-2.0.1')
    parser.add_argument('--classified', metavar='classified',
                        help='e.g., 1')

    args = parser.parse_args()

    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    sample_name = args.sample_name
    blast_program_dir = args.blast_program_dir
    RepeatModeler_Home = args.RepeatModeler_Home
    classified = args.classified

    log = Logger('HiTE.log', level='debug')

    tools_dir = os.getcwd() + '/../tools'

    confident_TE_path = tmp_output_dir + '/confident_TE.fa'
    confident_TE_consensus = tmp_output_dir + '/confident_TE.cons.fa'

    rename_fasta(confident_TE_path, confident_TE_path)
    contignames, contigs = read_fasta(confident_TE_path)
    new_contigs = {}
    for name in contignames:
        seq = contigs[name]
        if len(seq) < 100:
            continue
        new_contigs[name] = seq
    store_fasta(new_contigs, confident_TE_path)

    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_TE_path + ' -o ' + confident_TE_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    if classified is not None and int(classified) == 1:
        # # use RepeatClassifier to classify TE models
        # command = 'cd ' + tmp_output_dir + ' && ' + RepeatModeler_Home + '/RepeatClassifier -pa ' + str(threads) + ' -consensi ' + confident_TE_consensus
        # log.logger.debug(command)
        # os.system(command + ' > /dev/null 2>&1')

        TEClass_home = os.getcwd() + '/../classification'
        TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                          + ' --consensus ' + confident_TE_consensus + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(48) + ' -o ' + tmp_output_dir \
                          + ' --RepeatModeler_Home ' + RepeatModeler_Home
        log.logger.debug(TEClass_command)
        os.system(TEClass_command)

        # 把unknown的序列放在最后
        classified_TE_path = confident_TE_consensus + '.classified'
        names, contigs = read_fasta(classified_TE_path)
        names.sort(key=lambda x: x.split('#')[1])
        with open(classified_TE_path, 'w') as f_save:
            for name in names:
                f_save.write('>'+name+'\n'+contigs[name]+'\n')

