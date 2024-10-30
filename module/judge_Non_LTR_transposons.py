#-- coding: UTF-8 --
import argparse
import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, rename_fasta, Logger, file_exist, \
    flank_region_align_v5, get_candidate_non_ltr_parallel, get_domain_info

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE de novo Non_LTR module...')
    parser.add_argument('--seqs', metavar='seqs',
                        help='Please enter the result of de novo TE searching in HiTE, typically named longest_repeats_*.fa. Please provide the absolute path.')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--subset_script_path', metavar='subset_script_path',
                        help='Script to obtain a subset of high-copy TEs. Please provide the absolute path.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--library_dir', metavar='library_dir',
                        help='Please enter the directory for non_ltr library. Use an absolute path.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('--plant', metavar='plant',
                        help='Is it a plant genome, 1: true, 0: false.')
    parser.add_argument('--debug', metavar='recover',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false.')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='The flanking length of candidates to find the true boundaries.')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='The current split genome index.')
    parser.add_argument('--is_denovo_nonltr', metavar='is_denovo_nonltr',
                        help='Whether to detect non-ltr de novo, 1: true, 0: false.')
    parser.add_argument('-r', metavar='Reference path',
                        help='Input Reference path.')
    parser.add_argument('--prev_TE', metavar='prev_TE',
                        help='TEs fasta file that has already been identified. Please use the absolute path.')


    args = parser.parse_args()
    longest_repeats_flanked_path = args.seqs
    threads = int(args.t)
    subset_script_path = args.subset_script_path
    tmp_output_dir = args.tmp_output_dir
    library_dir = args.library_dir
    recover = args.recover
    plant = int(args.plant)
    debug = args.debug
    flanking_len = int(args.flanking_len)
    ref_index = args.ref_index
    is_denovo_nonltr = int(args.is_denovo_nonltr)
    reference = args.r
    prev_TE = args.prev_TE

    reference = os.path.realpath(reference)

    if debug is None:
        debug = 0
    else:
        debug = int(debug)

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()

    tmp_output_dir = os.path.abspath(tmp_output_dir)
    work_dir = tmp_output_dir

    log = Logger(tmp_output_dir + '/HiTE_Non_LTR.log', level='debug')

    TE_type = 'non_ltr'

    candidate_non_ltr_path = work_dir + '/candidate_non_ltr_' + str(ref_index) + '.fa'
    confident_non_ltr_path = work_dir + '/confident_non_ltr_' + str(ref_index) + '.fa'
    resut_file = confident_non_ltr_path
    if is_denovo_nonltr == 1:
        if not is_recover or not file_exist(resut_file):
            split_ref_dir = work_dir + '/ref_chr'
            LINE_domain_path = library_dir + '/LINEPeps.lib'
            remain_candidate_LINE_path = work_dir + '/remain_candidate_LINE_' + str(ref_index) + '.fa'
            remain_confident_LINE_path = work_dir + '/remain_confident_LINE_' + str(ref_index) + '.fa'
            output_table = remain_candidate_LINE_path + '.domain'
            temp_dir = work_dir + '/candidate_line_domain_' + str(ref_index)

            # Identify Non-LTR elements
            # 1. Identify candidate sequences with polyA/T + TSD structure from HiTE-FMEA results.
            candidate_SINE_path, candidate_LINE_path = get_candidate_non_ltr_parallel(longest_repeats_flanked_path, work_dir, threads)
            os.system('cat ' + candidate_SINE_path + ' > ' + candidate_non_ltr_path)
            os.system('cat ' + candidate_LINE_path + ' >> ' + candidate_non_ltr_path)
            # 2. Conduct homology search on candidate sequences, and search for polyA tails near homologous boundaries.
            flank_region_align_v5(candidate_non_ltr_path, confident_non_ltr_path, flanking_len, reference, split_ref_dir,
                                  TE_type, work_dir, threads, ref_index, log, subset_script_path,
                                  plant, debug, 0, result_type='cons')

            # 3. Select unrestrained LINE elements from candidate_LINE_path, align them to the domain, and extract reliable LINE elements.
            line_names, line_contigs = read_fasta(candidate_LINE_path)
            non_ltr_names, non_ltr_contigs = read_fasta(confident_non_ltr_path)
            remain_line_contigs = {}
            for name in line_names:
                if not non_ltr_contigs.__contains__(name):
                    remain_line_contigs[name] = line_contigs[name]
            store_fasta(remain_line_contigs, remain_candidate_LINE_path)
            get_domain_info(remain_candidate_LINE_path, LINE_domain_path, output_table, threads, temp_dir)
            domain_names, domain_contigs = read_fasta(LINE_domain_path)
            confident_LINE_contigs = {}
            name_set = set()
            with open(output_table, 'r') as f_r:
                for i, line in enumerate(f_r):
                    if i < 2:
                        continue
                    parts = line.split('\t')
                    # 如果包含完整的domain元素，则认为是真LINE元素
                    domain_name = parts[1]
                    domain_start = int(parts[4])
                    domain_end = int(parts[5])
                    if abs(domain_end-domain_start)/len(domain_contigs[domain_name]) >= 0.95:
                        name_set.add(parts[0])
            for name in name_set:
                confident_LINE_contigs[name] = line_contigs[name]
            store_fasta(confident_LINE_contigs, remain_confident_LINE_path)
            # 4. Add reliable LINE elements to confident_non_ltr_path and generate a consensus sequence.
            os.system('cat ' + remain_confident_LINE_path + ' >> ' + confident_non_ltr_path)

            confident_non_ltr_cons = confident_non_ltr_path + '.cons'
            cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                             + ' -G 0 -g 1 -A 80 -i ' + confident_non_ltr_path + ' -o ' + confident_non_ltr_cons + ' -T 0 -M 0' + ' > /dev/null 2>&1'
            os.system(cd_hit_command)
            rename_fasta(confident_non_ltr_cons, confident_non_ltr_path, 'Non-LTR_' + str(ref_index))
        else:
            log.logger.info(resut_file + ' exists, skip...')
    else:
        os.system('touch ' + confident_non_ltr_path)

    os.system('cat ' + resut_file + ' >> ' + prev_TE)