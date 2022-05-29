import codecs
import json
import os
import sys
import time

import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output, \
    get_multiple_alignment_repeat, split2cluster, cut_repeat_v1, judgeReduceThreads, get_ltr_suppl_from_ltrfinder, \
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v3, compare_seq, \
    getRegionCombination, getCombineFragments, convertToUpperCase_v1, generate_candidate_repeats_v2

def test(num):
    for i in range(10000):
        num += num
    return num

def cut_repeat(sam_path_bwa, repeats_path, cut_repeats_path):
    query_records = {}
    samfile = pysam.AlignmentFile(sam_path_bwa, "rb")
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        query_name = read.query_name
        reference_name = read.reference_name
        cigar = read.cigartuples
        cigarstr = read.cigarstring
        NM_tag = 0
        try:
            NM_tag = read.get_tag('NM')
        except KeyError:
            NM_tag = -1
        identity = compute_identity(cigarstr, NM_tag, 'BLAST')
        identity = float(identity) * 100
        is_reverse = read.is_reverse
        alignment_len = read.query_alignment_length
        # pos start from 1, change to 0
        q_start = int(read.query_alignment_start)  # [q_start, q_end)
        q_end = int(read.query_alignment_end)
        if q_start > q_end:
            tmp = q_start
            q_start = q_end
            q_end = tmp
        if not query_records.__contains__(query_name):
            query_records[query_name] = []
        records = query_records[query_name]
        records.append((reference_name, alignment_len, identity, q_start, q_end))
        query_records[query_name] = records

    repeat_contignames, repeat_contigs = read_fasta(repeats_path)
    cut_repeats = {}
    for query_name in query_records.keys():
        query_seq = repeat_contigs[query_name]
        query_len = len(query_seq)
        records = query_records[query_name]
        for i, record in enumerate(records):
            # filter first alignment
            if i == 0:
                continue
            identity = record[2]
            q_start = record[3]
            q_end = record[4]
            if identity < 95:
                continue
            # get repeats boundary by getting all alignment sequences
            new_seq = query_seq[q_start: q_end]
            new_query_name = query_name + '-p_' + str(i) + '-len_' + str(len(new_seq))
            cut_repeats[new_query_name] = new_seq
    store_fasta(cut_repeats, cut_repeats_path)

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
    # param_config_path = os.getcwd() + "/ParamConfig.json"
    # # read param config
    # with open(param_config_path, 'r') as load_f:
    #     param = json.load(load_f)
    # load_f.close()

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-26.10-28-8'
    reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'

    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)
    output_dir = tmp_output_dir
    repeats_consensus = tmp_output_dir + '/repeats.fa'
    unique_kmer_path = tmp_output_dir + '/kmer.txt'
    skip_threshold = 200
    identity_threshold = 0.90
    length_similarity_cutoff = 0.90
    tandem_region_cutoff = 0.5
    k_num = 31
    threads = 48
    partitions_num = threads
    tools_dir = os.getcwd() + '/tools'
    alias = 'dmel'
    # chrom_seg_length = int(param['chrom_seg_length'])
    # fault_tolerant_bases = 100

    reference = '/public/home/hpc194701009/KmerRepFinder_test/genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna'
    repbase = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepBase26.05.fasta/drorep.ref'
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(repbase, reference, use_align_tools, threads, tools_dir)

#
#
#     # Step1: generate candidate repeat regions
#     unique_kmer_map = {}
#     with open(unique_kmer_path, 'r') as f_r:
#         for line in f_r:
#             line = line.replace('\n', '')
#             kmer = line.split(' ')[0]
#             r_kmer = getReverseSequence(kmer)
#             unique_key = kmer if kmer < r_kmer else r_kmer
#             if unique_key.__contains__('N'):
#                 continue
#             unique_kmer_map[unique_key] = 1
#
#     # using multiple threads to gain speed
#     reference_pre = convertToUpperCase_v1(reference)
#     reference_tmp = multi_line(reference_pre, chrom_seg_length, k_num)
#
#     segments = []
#     with open(reference_tmp, 'r') as f_r:
#         for line in f_r:
#             line = line.replace('\n', '')
#             segments.append(line)
#     segments_cluster = split2cluster(segments, partitions_num)
#
#     ex = ProcessPoolExecutor(partitions_num)
#     repeat_dict = {}
#     jobs = []
#     for partiton_index in segments_cluster.keys():
#         cur_segments = segments_cluster[partiton_index]
#         job = ex.submit(generate_candidate_repeats_v2, cur_segments, k_num, unique_kmer_map, partiton_index,
#                         fault_tolerant_bases)
#         jobs.append(job)
#     ex.shutdown(wait=True)
#
#     for job in as_completed(jobs):
#         cur_repeat_dict = job.result()
#         for ref_name in cur_repeat_dict.keys():
#             parts = ref_name.split('$')
#             true_ref_name = parts[0]
#             start_pos = int(parts[1])
#             if not repeat_dict.__contains__(true_ref_name):
#                 repeat_dict[true_ref_name] = []
#             new_repeat_list = repeat_dict[true_ref_name]
#             cur_repeat_list = cur_repeat_dict[ref_name]
#             for repeat_item in cur_repeat_list:
#                 new_repeat_item = (start_pos + repeat_item[0], start_pos + repeat_item[1], repeat_item[2])
#                 new_repeat_list.append(new_repeat_item)
#     for ref_name in repeat_dict.keys():
#         repeat_list = repeat_dict[ref_name]
#         repeat_list.sort(key=lambda x: (x[1], x[2]))
#
#     repeats_path = tmp_output_dir + '/repeats.fa'
#     node_index = 0
#     with open(repeats_path, 'w') as f_save:
#         for ref_name in repeat_dict.keys():
#             repeat_list = repeat_dict[ref_name]
#             for repeat_item in repeat_list:
#                 start_pos = repeat_item[0]
#                 end_pos = repeat_item[1]
#                 query_name = 'N' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
#                 repeat = repeat_item[2]
#                 f_save.write('>' + query_name + '\n' + repeat + '\n')
#                 node_index += 1
#
#      # store repeat_dict for testing
#     repeat_dict_file = tmp_output_dir + '/repeat_dict.csv'
#     with codecs.open(repeat_dict_file, 'w', encoding='utf-8') as f:
#         json.dump(repeat_dict, f)
#
#     # Step2: ensure repeats boundary
#     # round-1
#     repeats_path = tmp_output_dir + '/repeats.fa'
#     use_align_tools = 'bwa'
#     sam_path_bwa = run_alignment(repeats_path, reference, use_align_tools, threads, tools_dir)
# #sam_path_bwa = tmp_output_dir + '/repeats.sam'
#     cut_repeats_path = tmp_output_dir + '/repeats.cut.fa'
#     cut_repeat(sam_path_bwa, repeats_path, cut_repeats_path)
#
#     # Step3: merge redundant sequences
#     cut_repeats_consensus = tmp_output_dir + '/repeats.cut.cons.fa'
#     cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.95) + ' -i ' + cut_repeats_path + ' -o ' + cut_repeats_consensus + ' -T 0 -M 0'
# #log.logger.debug(cd_hit_command)
#     os.system(cd_hit_command)
#
#
#     # Step4: get multiple alignment sequences
#     cut_repeats_consensus = tmp_output_dir + '/repeats.cut.cons.fa'
#     use_align_tools = 'bwa'
#     sam_path_bwa = run_alignment(cut_repeats_consensus, reference, use_align_tools, threads, tools_dir)
# #sam_path_bwa = tmp_output_dir + '/repeats.cut.cons.sam'
#     sam_paths = []
#     sam_paths.append(sam_path_bwa)
#     new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, cut_repeats_consensus)
#     repeat_multiple_path = tmp_output_dir + '/repeats.cut.cons.multiple.fa'
#     cut_repeat_contigNames, cut_repeat_contigs = read_fasta(cut_repeats_consensus)
#     node_index = 0
#     with open(repeat_multiple_path, 'w') as f_save:
#         for repeat_id in new_mapping_repeatIds.keys():
#             freq = new_mapping_repeatIds[repeat_id][0]
#             seq = cut_repeat_contigs[repeat_id]
#             #f_save.write('>' + repeat_id + '\n' + seq + '\n')
#             #f_save.write('>' + repeat_id + '\tcopies=' + str(freq) + '\n' + seq + '\n')
#             f_save.write('>N' + str(node_index) + '\n' + seq + '\n')
#             node_index += 1

#     repeat_multiple_path = tmp_output_dir + '/repeats.cut.cons.multiple.fa'
#     starttime = time.time()
#     # Step0. use RepeatMasker/trf to mask all low complexity/tandem repeats in raw repeat region
#     # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
#     RepeatMasker_Home = param['RepeatMasker_Home']
#     RepeatMasker_output_dir = tmp_output_dir + '/noint'
#     RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
#                            + ' -noint -x -dir ' + RepeatMasker_output_dir + ' ' + repeat_multiple_path
# #os.system('rm -rf ' + RepeatMasker_output_dir)
# #log.logger.debug(RepeatMasker_command)
#     #os.system(RepeatMasker_command)
#
#     (repeat_multiple_dir, repeat_multiple_filename) = os.path.split(repeat_multiple_path)
#     (repeat_multiple_name, repeat_multiple_extension) = os.path.splitext(repeat_multiple_filename)
#     trf_masked_repeats = RepeatMasker_output_dir + '/' + repeat_multiple_filename + '.masked'
#     trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
#     repeats_contigNames, repeats_contigs = read_fasta(repeat_multiple_path)
#     repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
#     with open(repeats_path, 'w') as f_save:
#         for name in trf_contigNames:
#             seq = trf_contigs[name]
#             if float(seq.count('X')) / len(seq) < tandem_region_cutoff:
#                 f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
#     RepeatMasker_Home = param['RepeatMasker_Home']
#     RepeatMasker_output_dir = tmp_output_dir + '/noint'
#     RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(
#         threads) + ' -lib ' + repeats_path + ' -dir ' + RepeatMasker_output_dir + ' ' + reference
#     os.system(RepeatMasker_command)

#     TRF_Path = param['TRF_Path']
#
#     trf_dir = tmp_output_dir + '/trf_temp'
#     if not os.path.exists(trf_dir):
#         os.makedirs(trf_dir)
#
#     trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + reference + ' 2 7 7 80 10 50 500 -f -d -m'
# #log.logger.debug(trf_command)
#     os.system(trf_command)
#     trf_masked_repeats = trf_dir + '/' + ref_filename + '.2.7.7.80.10.50.500.mask'

    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    # repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # with open(repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')

    # endtime = time.time()
    # dtime = endtime - starttime
#log.logger.debug("Step0: use trf to mask genome: %.8s s" % (dtime))











    #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/krf_output/CRD.2022-05-25.20-12-12'

    # Step3: According to valid paths, the fragments in the region are connected across gaps.
    # Fragments in region are connected according to the longer path in valid paths.


    # connected_frags_file = tmp_output_dir + '/connected_frags.csv'
    # file = open(connected_frags_file, 'r')
    # js = file.read()
    # connected_frags = json.loads(js)
    #
    # refNames, refContigs = read_fasta(reference)
    # repeats_connected_file = tmp_output_dir + '/repeats_connected.fa'
    # repeats_connected = {}
    # index = 0
    # for region_index in connected_frags.keys():
    #     for connected_frag in connected_frags[region_index]:
    #         frag_name = connected_frag[0].split(',')[0]
    #         ref_name = frag_name.split('-s_')[1].split('-')[0]
    #         query_name = 'R' + str(index) + '-' + frag_name
    #         seq = refContigs[ref_name][connected_frag[1]: connected_frag[2] + 1]
    #         index += 1
    #         repeats_connected[query_name] = seq
    # sorted_repeats_connected = {k: v for k, v in sorted(repeats_connected.items(), key=lambda item: -len(item[1]))}
    # store_fasta(sorted_repeats_connected, repeats_connected_file)