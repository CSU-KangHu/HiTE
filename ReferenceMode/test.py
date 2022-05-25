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

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-24.9-36-59'
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

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/krf_output/CRD.2022-05-25.20-12-12'

    # Step3: According to valid paths, the fragments in the region are connected across gaps.
    # Fragments in region are connected according to the longer path in valid paths.

    region_paths_file = tmp_output_dir + '/region_paths.csv'
    file = open(region_paths_file, 'r')
    js = file.read()
    region_paths = json.loads(js)
    region_fragments_file = tmp_output_dir + '/region_fragments.csv'
    file = open(region_fragments_file, 'r')
    js = file.read()
    region_fragments = json.loads(js)

    # keeped_path = {region_id: [f1f2f3, f4, f5f6]}
    keeped_paths = {}
    for region_index in region_paths.keys():
        region_fragment = region_fragments[region_index]
        for path in region_paths[region_index]:
            isPathExist = True
            pathName = ''
            for i, frag_index in enumerate(path.keys()):
                if region_fragment.__contains__(frag_index):
                    frag_item = region_fragment[frag_index]
                    frag_name = frag_item[0]
                    if i == 0:
                        pathName += frag_name
                    else:
                        pathName += ',' + frag_name
                else:
                    isPathExist = False
                    break

            if isPathExist:
                if not keeped_paths.__contains__(region_index):
                    keeped_paths[region_index] = []
                paths = keeped_paths[region_index]
                paths.append(pathName)
                for frag_index in path.keys():
                    del region_fragment[frag_index]

    # store keeped_paths for testing
    keeped_paths_file = tmp_output_dir + '/keeped_paths.csv'
    with codecs.open(keeped_paths_file, 'w', encoding='utf-8') as f:
        json.dump(keeped_paths, f)