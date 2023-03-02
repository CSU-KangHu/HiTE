import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta_v1, store_fasta

def get_perfect_name(perfect_path):
    perfect_names = set()
    with open(perfect_path, 'r') as f_r:
        for line in f_r:
            name = line.replace('\n', '')
            perfect_names.add(name)
    f_r.close()
    return perfect_names

if __name__ == '__main__':
    #我们获得的valid_tir.fa.itr.filter_dup只能获得335个perfect families，而valid_tir.fa(直接使用最长的lenTE)得到481个perfect families。我们需要找一下哪些是我们丢失的。
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library'

    test_dir_273 = tmp_dir + '/get_family_summary_test'
    perfect_273_path = test_dir_273 + '/perfect.families'
    test_432_dir = tmp_dir + '/get_family_summary_test_522'
    perfect_432_path = test_432_dir + '/perfect.families'

    perfect_432 = get_perfect_name(perfect_432_path)
    perfect_273 = get_perfect_name(perfect_273_path)

    diff_set = perfect_432.difference(perfect_273)
    print(diff_set)
    print(len(diff_set))
    #335 perfect.families