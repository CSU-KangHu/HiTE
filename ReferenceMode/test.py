import codecs
import json
import os
import sys
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output, \
    get_multiple_alignment_repeat, split2cluster, cut_repeat_v1, judgeReduceThreads, get_ltr_suppl_from_ltrfinder, \
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v3, compare_seq, getRegionCombination, getCombineFragments

def test(num):
    for i in range(10000):
        num += num
    return num

if __name__ == '__main__':
    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-21.0-9-18/region_combination_tmp'
    # file = tmp_output_dir + '/rc_0.csv'
    # region_combination_size = float(sys.getsizeof(region_combination))

    # ex = ProcessPoolExecutor(48)
    # bigint = 1024*1024*1024
    # big_list = list(x for x in range(bigint))
    #
    # MAX_JOBS_IN_QUEUE = 500
    # jobs_left = len(big_list)
    # jobs_iter = iter(big_list)
    # jobs = {}
    # while jobs_left:
    #     for num in jobs_iter:
    #         job = ex.submit(test, num)
    #         jobs[job] = 1
    #         if len(jobs) > MAX_JOBS_IN_QUEUE:
    #             break  # limit the job submission for now job
    #
    #     for job in as_completed(jobs):
    #         jobs_left -= 1
    #         num = job.result()
    #         del jobs[job]
    #         break
    # ex.shutdown(wait=True)
    param_config_path = os.getcwd() + "/ParamConfig.json"
    # read param config
    with open(param_config_path, 'r') as load_f:
        param = json.load(load_f)
    load_f.close()

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-26.10-28-8'
    reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'
    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)
    output_dir = tmp_output_dir
    repeats_consensus = tmp_output_dir + '/repeats.fa'
    skip_threshold = 200
    identity_threshold = 0.90
    length_similarity_cutoff = 0.90
    tandem_region_cutoff = 0.8
    k_num = 31*3
    threads = 48
    partitions_num = threads
    tools_dir = os.getcwd() + '/tools'
    alias = 'dmel'

    #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/krf_output/CRD.2022-05-25.20-12-12'

    # Step3: According to valid paths, the fragments in the region are connected across gaps.
    # Fragments in region are connected according to the longer path in valid paths.


    connected_frags_file = tmp_output_dir + '/connected_frags.csv'
    file = open(connected_frags_file, 'r')
    js = file.read()
    connected_frags = json.loads(js)

    refNames, refContigs = read_fasta(reference)
    repeats_connected_file = tmp_output_dir + '/repeats_connected.fa'
    repeats_connected = {}
    index = 0
    for region_index in connected_frags.keys():
        for connected_frag in connected_frags[region_index]:
            frag_name = connected_frag[0].split(',')[0]
            ref_name = frag_name.split('-s_')[1].split('-')[0]
            query_name = 'R' + str(index) + '-' + frag_name
            seq = refContigs[ref_name][connected_frag[1]: connected_frag[2] + 1]
            index += 1
            repeats_connected[query_name] = seq
    sorted_repeats_connected = {k: v for k, v in sorted(repeats_connected.items(), key=lambda item: -len(item[1]))}
    store_fasta(sorted_repeats_connected, repeats_connected_file)