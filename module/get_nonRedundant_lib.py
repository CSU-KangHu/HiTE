import argparse
import os
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, rename_fasta, remove_ltr_from_tir



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
    parser.add_argument('--test_home', metavar='test_home',
                        help='The directory of HiTE module')

    args = parser.parse_args()

    threads = int(args.t)
    confident_ltr_cut_path = args.confident_ltr_cut
    confident_tir_path = args.confident_tir
    confident_helitron_path = args.confident_helitron
    confident_non_ltr_path = args.confident_non_ltr
    confident_other_path = args.confident_other
    tmp_output_dir = args.tmp_output_dir
    test_home = args.test_home

    confident_ltr_cut_path = os.path.realpath(confident_ltr_cut_path)
    confident_tir_path = os.path.realpath(confident_tir_path)
    confident_helitron_path = os.path.realpath(confident_helitron_path)
    confident_non_ltr_path = os.path.realpath(confident_non_ltr_path)
    confident_other_path = os.path.realpath(confident_other_path)

    tmp_output_dir = os.path.abspath(tmp_output_dir)

    log = Logger(tmp_output_dir+'/HiTE_lib.log', level='debug')

    log.logger.info(confident_ltr_cut_path)
    log.logger.info(confident_tir_path)
    log.logger.info(confident_helitron_path)
    log.logger.info(confident_non_ltr_path)
    log.logger.info(confident_other_path)
    log.logger.info(tmp_output_dir)

    rename_fasta(confident_tir_path, confident_tir_path, 'TIR')
    rename_fasta(confident_helitron_path, confident_helitron_path, 'Helitron')
    rename_fasta(confident_non_ltr_path, confident_non_ltr_path, 'Denovo_Non_LTR')


    final_confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    final_confident_helitron_path = tmp_output_dir + '/confident_helitron.fa'
    final_confident_non_ltr_path = tmp_output_dir + '/confident_non_ltr.fa'
    final_confident_other_path = tmp_output_dir + '/confident_other.fa'

    # clustering and rename
    candidate_tir_cons = confident_tir_path + '.cons'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_tir_path + ' -o ' + candidate_tir_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)
    rename_fasta(candidate_tir_cons, final_confident_tir_path, 'TIR')

    confident_helitron_cons = confident_helitron_path + '.cons'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_helitron_path + ' -o ' + confident_helitron_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)
    rename_fasta(confident_helitron_cons, final_confident_helitron_path, 'Helitron')

    confident_non_ltr_cons = confident_non_ltr_path + '.cons'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_non_ltr_path + ' -o ' + confident_non_ltr_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)
    rename_fasta(confident_non_ltr_cons, final_confident_non_ltr_path, 'Denovo_Non_LTR')

    # Remove TIR elements containing LTR
    remove_ltr_from_tir(confident_ltr_cut_path, final_confident_tir_path, threads)

    # Merge all TE types (TIR+Helitron+Non_LTR+Other)
    confident_TE_path = tmp_output_dir + '/confident_TE.fa'
    os.system('cat ' + final_confident_tir_path + ' > ' + confident_TE_path)
    os.system('cat ' + final_confident_helitron_path + ' >> ' + confident_TE_path)
    os.system('cat ' + final_confident_non_ltr_path + ' >> ' + confident_TE_path)
    os.system('cat ' + confident_other_path + ' >> ' + confident_TE_path)
    os.system('cat ' + confident_ltr_cut_path + ' >> ' + confident_TE_path)

    # Unpack nested TEs within TEs
    clean_TE_path = tmp_output_dir + '/confident_TE.clean.fa'
    remove_nested_command = 'python3 ' + test_home + '/remove_nested_lib.py ' \
                            + ' -t ' + str(threads) \
                            + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(5) \
                            + ' --input1 ' + confident_TE_path \
                            + ' --input2 ' + confident_TE_path \
                            + ' --output ' + clean_TE_path
    os.system(remove_nested_command)

    sample_name = 'test'
    confident_TE_consensus = tmp_output_dir + '/confident_TE.cons.fa'

    #rename_fasta(clean_TE_path, clean_TE_path)
    contignames, contigs = read_fasta(clean_TE_path)
    new_contigs = {}
    for name in contignames:
        seq = contigs[name]
        if len(seq) < 100:
            continue
        new_contigs[name] = seq
    store_fasta(new_contigs, clean_TE_path)

    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + clean_TE_path + ' -o ' + confident_TE_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)



