#-- coding: UTF-8 --
import argparse
import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, getReverseSequence, read_fasta_v2, \
    multi_process_align_and_get_copies, rename_fasta, Logger, file_exist


def extract_sequence_from_db(db_path, features, store_path):
    store_contigs = {}
    db_contignames, db_contigs = read_fasta(db_path)
    for name in db_contignames:
        for feature in features:
            if name.__contains__(feature):
                db_seq = db_contigs[name]
                store_contigs[name] = db_seq
    store_fasta(store_contigs, store_path)

def preprocess():
    db_path = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker/Libraries/RepeatMasker.lib'
    features = ['#LINE', '#SINE', 'DIRS', '#DNA/Crypton']
    store_path = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/non_LTR.lib'
    extract_sequence_from_db(db_path, features, store_path)

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE homology Non_LTR module...')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--library_dir', metavar='library_dir',
                        help='Please enter the directory for non_ltr library. Use an absolute path.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('-r', metavar='Reference path',
                        help='Input Reference path')


    args = parser.parse_args()
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    library_dir = args.library_dir
    recover = args.recover
    reference = args.r

    reference = os.path.realpath(reference)

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()

    tmp_output_dir = os.path.abspath(tmp_output_dir)

    log = Logger(tmp_output_dir + '/HiTE_other.log', level='debug')

    confident_other_path = tmp_output_dir + '/confident_other.fa'
    resut_file = confident_other_path
    if not is_recover or not file_exist(resut_file):
        non_LTR_lib = library_dir + '/non_LTR.lib'
        other_TE_dir = tmp_output_dir + '/other_TE'
        os.system('rm -rf ' + other_TE_dir)
        if not os.path.exists(other_TE_dir):
            os.makedirs(other_TE_dir)

        confident_non_ltr_contigs = {}
        # 1. Align the library to the reference to obtain copies.
        TE_type = 'non_ltr'
        all_copies = multi_process_align_and_get_copies(non_LTR_lib, reference, other_TE_dir, TE_type, threads, is_removed_dir=True, query_coverage=0.95, subject_coverage=0)

        # 2. Take the longest copy as the identified non-LTR element.
        ref_names, ref_contigs = read_fasta(reference)
        new_all_copies = {}
        for query_name in all_copies.keys():
            copies = all_copies[query_name]
            max_len_copy_name = ''
            max_len = 0
            for copy in copies:
                ref_name = copy[0]
                copy_ref_start = int(copy[1])
                copy_ref_end = int(copy[2])
                copy_len = copy_ref_end - copy_ref_start + 1
                if copy_ref_start - 1 < 0 or copy_ref_end > len(ref_contigs[ref_name]):
                    continue
                copy_seq = ref_contigs[ref_name][copy_ref_start - 1: copy_ref_end]
                if len(copy_seq) < 100:
                    continue
                new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end)
                if len(copy_seq) > max_len:
                    confident_non_ltr_contigs[new_name] = copy_seq
                    max_len = len(copy_seq)
        store_fasta(confident_non_ltr_contigs, confident_other_path)
        rename_fasta(confident_other_path, confident_other_path, 'Homology_Non_LTR')
    else:
        log.logger.info(resut_file + ' exists, skip...')
