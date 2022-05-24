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

    # new strategy by Kang Hu 2022/05/24
    # Step1: generate repeats.fa and connected_regions
    # load repeat_dict
    repeat_dict_file = tmp_output_dir + '/repeat_dict.csv'
    file = open(repeat_dict_file, 'r')
    js = file.read()
    repeat_dict = json.loads(js)

    # connected_regions = {ref_name: {region_id: [(f1, start1, end1), (f2, start2, end2), (f3, start3, end3)], [(f4, start4, end4), (f5, start5, end5), (f6, start6, end6)]}}
    connected_regions = {}
    repeats_path = tmp_output_dir + '/repeats.fa'
    node_index = 0
    region_index = 0
    with open(repeats_path, 'w') as f_save:
        for ref_name in repeat_dict.keys():
            repeat_list = repeat_dict[ref_name]
            if not connected_regions.__contains__(ref_name):
                connected_regions[ref_name] = {}
            regions = connected_regions[ref_name]
            last_start_pos = -1
            last_end_pos = -1
            for repeat_item in repeat_list:
                start_pos = repeat_item[0]
                end_pos = repeat_item[1]
                query_name = 'N' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
                repeat = repeat_item[2]
                f_save.write('>' + query_name + '\n' + repeat + '\n')
                node_index += 1
                # generate connected_regions
                if last_start_pos == -1:
                    regions[region_index] = [(query_name, start_pos, end_pos)]
                else:
                    if (start_pos - last_end_pos) < skip_threshold:
                        # close to current region
                        cur_region = regions[region_index]
                        cur_region.append((query_name, start_pos, end_pos))
                        regions[region_index] = cur_region
                    else:
                        # far from current region, start a new region
                        region_index += 1
                        cur_region = []
                        cur_region.append((query_name, start_pos, end_pos))
                        regions[region_index] = cur_region
                last_start_pos = start_pos
                last_end_pos = end_pos
            connected_regions[ref_name] = regions

    # store connected_regions for testing
    connected_regions_file = tmp_output_dir + '/connected_regions.csv'
    with codecs.open(connected_regions_file, 'w', encoding='utf-8') as f:
        json.dump(connected_regions, f)

    # Step2: align repeat back to reference, generate frag_pos_dict
    repeat_contignames, repeat_contigs = read_fasta(repeats_path)
    use_align_tools = 'bwa'
    # sam_path_bwa = run_alignment(repeats_path, reference, use_align_tools, threads, tools_dir)
    sam_path_bwa = tmp_output_dir + '/repeats.sam'
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
        q_start = int(read.query_alignment_start)
        q_end = int(read.query_alignment_end)
        t_start = int(read.reference_start)
        t_end = int(read.reference_end)
        if t_start > t_end:
            tmp = t_start
            t_start = t_end
            t_end = tmp

        if not query_records.__contains__(query_name):
            query_records[query_name] = []
        records = query_records[query_name]
        records.append((reference_name, alignment_len, identity, t_start, t_end))
        query_records[query_name] = records

    # frag_pos_dict = {f1: {ref1: [(start1, end1), (start2, end2), (start3, end3)]}}
    frag_pos_dict = {}
    for query_name in query_records.keys():
        complete_alignment_num = 0
        query_seq = repeat_contigs[query_name]
        query_len = len(query_seq)
        if not frag_pos_dict.__contains__(query_name):
            frag_pos_dict[query_name] = {}
        ref_pos = frag_pos_dict[query_name]
        # get fragments original position
        pos_parts = query_name.split('-s_')[1].split('-')
        original_ref = pos_parts[0]
        original_start = int(pos_parts[1])
        original_end = int(pos_parts[2])
        for i, record in enumerate(query_records[query_name]):
            reference_name = record[0]
            alignment_len = record[1]
            identity = record[2]
            t_start = record[3]
            t_end = record[4]
            if float(alignment_len) / query_len >= 0.95 and identity >= 95:
                complete_alignment_num += 1
            if not ref_pos.__contains__(reference_name):
                ref_pos[reference_name] = []
            same_chr_pos = ref_pos[reference_name]
            # alignment from original parts
            if original_ref == reference_name and abs(t_start-original_start) < 10 and abs(t_end-original_end) < 10:
                continue
            same_chr_pos.append((t_start, t_end))
            ref_pos[reference_name] = same_chr_pos
        # sort pos in each ref_name
        for reference_name in ref_pos.keys():
            same_chr_pos = ref_pos[reference_name]
            same_chr_pos.sort(key=lambda x: (x[0], x[1]))
            ref_pos[reference_name] = same_chr_pos
        frag_pos_dict[query_name] = ref_pos

    # store frag_pos_dict for testing
    frag_pos_dict_file = tmp_output_dir + '/frag_pos_dict.csv'
    with codecs.open(frag_pos_dict_file, 'w', encoding='utf-8') as f:
        json.dump(frag_pos_dict, f)

    # Step3: generate pathMatrix
    # pathMatrix = {region_id: {ref_name: [[(start1, end1), (start2, end2), (start3, end3)], [(start1, end1), (start2, end2), (start3, end3)]]}}
    pathMatrix = {}
    for ref_name in connected_regions.keys():
        regions = connected_regions[ref_name]
        for region_index in regions.keys():
            if not pathMatrix.__contains__(region_index):
                pathMatrix[region_index] = {}
            cur_region_matrixs = pathMatrix[region_index]

            # get one region
            cur_region = regions[region_index]
            # get all ref_names for one region
            cur_ref_names_union = set()
            for i in range(len(cur_region)):
                frag_name = cur_region[i][0]
                if not frag_pos_dict.__contains__(frag_name):
                    frag_pos_dict[frag_name] = {}
                ref_pos = frag_pos_dict[frag_name]
                for ref in ref_pos.keys():
                    cur_ref_names_union.add(ref)

            for ref in cur_ref_names_union:
                if not cur_region_matrixs.__contains__(ref):
                    cur_region_matrixs[ref] = []
                cur_ref_matrix = cur_region_matrixs[ref]
                for i in range(len(cur_region)):
                    frag_name = cur_region[i][0]
                    if not frag_pos_dict.__contains__(frag_name):
                        frag_pos_dict[frag_name] = {}
                    ref_pos = frag_pos_dict[frag_name]
                    if not ref_pos.__contains__(ref):
                        ref_pos[ref] = []
                    same_chr_pos = ref_pos[ref]
                    cur_ref_matrix.append(same_chr_pos)
                cur_region_matrixs[ref] = cur_ref_matrix
            pathMatrix[region_index] = cur_region_matrixs


    # store pathMatrix for testing
    pathMatrix_file = tmp_output_dir + '/pathMatrix.csv'
    with codecs.open(pathMatrix_file, 'w', encoding='utf-8') as f:
        json.dump(pathMatrix, f)

    # load pathMatrix
    pathMatrix_file = tmp_output_dir + '/pathMatrix.csv'
    file = open(pathMatrix_file, 'r')
    js = file.read()
    pathMatrix = json.loads(js)

    # go through each Matrix, compute the longest path
    # specific_region_id = "1"
    # specific_ref_name = 'U'
    # for region_index in pathMatrix:
    #     if region_index == specific_region_id:
    #         cur_region_matrixs = pathMatrix[region_index]
    #         for ref in cur_region_matrixs.keys():
    #             if ref == specific_ref_name:
    #                 cur_ref_matrix = cur_region_matrixs[ref]
    #                 print(cur_ref_matrix)




    # unique_kmer_path = tmp_output_dir + '/kmer.txt'
    # freq_distribution = {}
    # total_kmer_num = 0
    # with open(unique_kmer_path, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         parts = line.split(' ')
    #         kmer = parts[0]
    #         r_kmer = getReverseSequence(kmer)
    #         unique_key = kmer if kmer < r_kmer else r_kmer
    #         freq = int(parts[1])
    #         if not freq_distribution.__contains__(freq):
    #             freq_distribution[freq] = 0
    #         num = freq_distribution[freq]
    #         freq_distribution[freq] = num + 1
    #         total_kmer_num += 1
    #
    # freq_distribution = {k: v for k, v in
    #                          sorted(freq_distribution.items(), key=lambda item: item[0])}
    # print(freq_distribution)
    # cur_kmer_num = 0
    # for repeat_num in freq_distribution:
    #     kmer_num = freq_distribution[repeat_num]
    #     cur_kmer_num += kmer_num
    #     if cur_kmer_num > total_kmer_num/2:
    #         print(repeat_num)
    #         break

    # # try to chain all fragments
    # # --------------------------------------------------------------------------------------
    # # newly strategy: 2022-04-29 by Kang Hu
    # # 01: use bwa to get single mapped sequence
    # candidate_repeats_path = repeats_consensus
    # blast_program_dir = param['RMBlast_Home']
    # use_align_tools = 'bwa'
    # sam_path_bwa = run_alignment(candidate_repeats_path, reference, use_align_tools, threads, tools_dir)
    # sam_paths = []
    # sam_paths.append(sam_path_bwa)
    # # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info(sam_paths)
    # new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, candidate_repeats_path)
    # repeat_freq_path = tmp_output_dir + '/repeats.freq.fa'
    # repeat_multiple_path = tmp_output_dir + '/repeats.multiple.fa'
    # merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # # single repeat probably be Chimeric repeat
    # with open(repeat_freq_path, 'w') as f_save:
    #     for repeat_id in new_mapping_repeatIds.keys():
    #         freq = new_mapping_repeatIds[repeat_id][0]
    #         seq = merge_repeat_contigs[repeat_id]
    #         if freq <= 1 or (freq < 5 and len(seq) < 80):
    #             continue
    #         f_save.write('>' + repeat_id + '\tcopies=' + str(freq) + '\n' + seq + '\n')
    #
    # with open(repeat_multiple_path, 'w') as f_save:
    #     for repeat_id in new_mapping_repeatIds.keys():
    #         freq = new_mapping_repeatIds[repeat_id][0]
    #         seq = merge_repeat_contigs[repeat_id]
    #         if freq <= 1:
    #             continue
    #         f_save.write('>' + repeat_id + '\n' + seq + '\n')
    #
    # # 06: merge
    # merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    # os.system('cat ' + repeat_multiple_path + ' > ' + merge_pure)
    # ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    # os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    # # cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -c ' + str(0.8) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # print(cd_hit_command)
    # os.system(cd_hit_command)
    #
    #
    # # Step0. use RepeatMasker/trf to mask all low complexity/tandem repeats in raw repeat region
    # # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    # TRF_Path = param['TRF_Path']
    #
    # trf_dir = tmp_output_dir + '/trf_temp'
    # if not os.path.exists(trf_dir):
    #     os.makedirs(trf_dir)
    # (repeat_dir, repeat_filename) = os.path.split(merge_pure_consensus)
    # (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    # trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + merge_pure_consensus + ' 2 7 7 80 10 50 500 -f -d -m'
    # print(trf_command)
    # os.system(trf_command)
    # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'
    #
    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    # repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # with open(repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
    #
    # # --------------------------------------------------------------------------------------
    # # Step10. run TE classification to classify TE family
    # print('Start step8: get classified consensus sequence')
    # sample_name = alias
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + repeats_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    # print(TEClass_command)
    # os.system(TEClass_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step11. assign a family name for each classified TE consensus
    # classified_consensus_path = repeats_path + '.final.classified'
    # classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    # family_path = tmp_output_dir + '/family_' + sample_name + '.fasta'
    # with open(family_path, 'w') as f_save:
    #     for f_id, name in enumerate(classified_contigNames):
    #         sequence = classified_contigs[name]
    #         class_name = name.split('#')[1]
    #         if len(sequence) < 80 and class_name == 'Unknown':
    #             continue
    #         f_save.write('>family-' + str(f_id) + '#' + class_name + '\n' + sequence + '\n')


    # reference = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/genome_model.fa'
    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/krf_output/CRD.2022-05-16.10-45-46'
    # merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    # tools_dir = os.getcwd() + '/tools'
    # cd_hit_command = tools_dir + '/cd-hit-est -s 0.95 -c 0.95 -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    # TRF_Path = '/public/home/hpc194701009/repeat_detect_tools/trf409.linux64'
    # tandem_region_cutoff = 0.9
    # trf_dir = tmp_output_dir + '/trf_temp'
    # if not os.path.exists(trf_dir):
    #     os.makedirs(trf_dir)
    # (repeat_dir, repeat_filename) = os.path.split(merge_pure_consensus)
    # (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    # trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + merge_pure_consensus + ' 2 7 7 80 10 50 500 -f -d -m'
    # os.system(trf_command)
    # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'
    #
    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    # repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # with open(repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
    #
    # # --------------------------------------------------------------------------------------
    # # Step10. run TE classification to classify TE family
    # sample_name = 'model'
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + repeats_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(48) + ' -o ' + tmp_output_dir
    # os.system(TEClass_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step11. assign a family name for each classified TE consensus
    # classified_consensus_path = repeats_path + '.final.classified'
    # classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    # family_path = tmp_output_dir + '/family_' + sample_name + '.fasta'
    # with open(family_path, 'w') as f_save:
    #     for f_id, name in enumerate(classified_contigNames):
    #         sequence = classified_contigs[name]
    #         # if len(sequence) < 80:
    #         #     continue
    #         class_name = name.split('#')[1]
    #         f_save.write('>family-' + str(f_id) + '#' + class_name + '\n' + sequence + '\n')
    #
    # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-16.17-18-37'
    # family_path = output_dir + '/family_dmel.fasta'
    # RepeatMasker_Home = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker'
    # RepeatMasker_command = 'cd ' + output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(48) \
    #                     + ' -noint -x ' + ' ' + reference
    # os.system('rm -rf ' + RepeatMasker_output_dir)
    # os.system(RepeatMasker_command)

    # curated_library = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/O.sativa_curated.fasta'
    # curated_contignames, curated_contigs = read_fasta(curated_library)
    # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/no_simple_repeats'
    # # # extract from curated library
    # model_lib_path = output_dir + '/O.sativa_curated.fasta'
    # with open(model_lib_path, 'w') as f_save:
    #     for name in curated_contignames:
    #         class_name = name.split('#')[1]
    #         if class_name == 'Simple_repeat' or class_name == 'Low_complexity':
    #             continue
    #         seq = curated_contigs[name]
    #         f_save.write('>' + name + '\n' + seq + '\n')

    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-16.17-18-37'
    # family_path = tmp_output_dir + '/family_dmel.fasta'
    # family_contigNames, family_contigs = read_fasta(family_path)
    # new_family_path = tmp_output_dir + '/family_dmel.filter_80.fasta'
    # with open(new_family_path, 'w') as f_save:
    #     for f_id, name in enumerate(family_contigNames):
    #         sequence = family_contigs[name]
    #         class_name = name.split('#')[1]
    #         if len(sequence) < 80:
    #             continue
    #         f_save.write('>family-' + str(f_id) + '#' + class_name + '\n' + sequence + '\n')

    # blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    # protein_db_path = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker/Libraries/RepeatPeps.lib'
    # # Step2. blastx build protein db
    # blast_db_command = blast_program_dir + '/bin/makeblastdb -dbtype prot -in ' + protein_db_path
    # print(blast_db_command)
    # #os.system(blast_db_command)
    #
    # families_path = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-13.10-18-52/repeats.filter_tandem.fa.final.classified'
    # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-13.10-18-52'
    # thread_num = '8'
    # # Step3. blastx align families to protein db
    # blastxResults_path = output_dir + '/tmpBlastXResults.out'
    # blast_align_command = blast_program_dir + '/bin/blastx -db ' + protein_db_path + ' -num_threads ' + thread_num + ' -query ' + families_path + ' -word_size 2 -outfmt 6 > ' + blastxResults_path
    # print(blast_align_command)
    # #os.system(blast_align_command)
    #
    # # parse blast output
    # high_confidence_results = {}
    # with open(blastxResults_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         family_id = parts[0]
    #         domain_id = parts[1]
    #         identity = float(parts[2])
    #         query_start = int(parts[6])
    #         query_end = int(parts[7])
    #         domain_start = int(parts[8])
    #         domain_end = int(parts[9])
    #         if identity >= 80:
    #             if not high_confidence_results.__contains__(family_id):
    #                 high_confidence_results[family_id] = []
    #             protein_alignments = high_confidence_results[family_id]
    #             protein_alignments.append((family_id, query_start, query_end, domain_id, domain_start, domain_end, identity))
    #             high_confidence_results[family_id] = protein_alignments
    #
    #
    # trust_TE = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-13.10-18-52/trust_TE.fa'
    # trust_TE_out = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-13.10-18-52/trust_TE_struct.out'
    # novel_TE = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-13.10-18-52/novel_TE.fa'
    #
    # contigNames, contigs = read_fasta(families_path)
    # with open(trust_TE, 'w') as f_save:
    #     for family_id in high_confidence_results.keys():
    #         f_save.write('>'+family_id+'\n'+contigs[family_id]+'\n')
    #
    # with open(novel_TE, 'w') as f_save:
    #     for family_id in contigNames:
    #         if family_id not in high_confidence_results.keys():
    #             f_save.write('>'+family_id+'\n'+contigs[family_id]+'\n')
    #
    # with open(trust_TE_out, 'w') as f_save:
    #     for family_id in high_confidence_results.keys():
    #         item_list = high_confidence_results[family_id]
    #         for item in item_list:
    #             f_save.write(str(item[0])+'\t'+str(item[1])+'\t'+str(item[2])+'\t'+str(item[3])+'\t'+str(item[4])+'\t'+str(item[5])+'\t'+str(item[6])+'\n')





    # #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/drerio/CRD.2022-05-10.2-32-16'
    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/oryza_sative/CRD.2022-05-11.11-45-41'
    # candidate_repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # candidate_consensus_path = tmp_output_dir + '/repeats.filter_tandem.consensus.fa'
    # candidate_contignames, candidate_contigs = read_fasta(candidate_consensus_path)
    #
    # new_candidate_contigs = {}
    # for name in candidate_contignames:
    #     seq = candidate_contigs[name]
    #     new_candidate_contigs[name.split('-')[0]] = seq
    #
    #
    # tools_dir = os.getcwd() + '/tools'
    # cd_hit_command = tools_dir + '/cd-hit-est -s 0.8 -c 0.8 -i ' + candidate_repeats_path + ' -o ' + candidate_consensus_path + ' -T 0 -M 0'
    # #os.system(cd_hit_command)
    # cluster_file = candidate_consensus_path + '.clstr'
    # cluster_info = {}
    # cluster_id = ''
    # with open(cluster_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         if line.startswith('>Cluster'):
    #             cluster_id = line
    #             continue
    #         if cluster_id != '':
    #             if not cluster_info.__contains__(cluster_id):
    #                 cluster_info[cluster_id] = []
    #             cluster_records = cluster_info[cluster_id]
    #             cluster_records.append(line)
    #             cluster_info[cluster_id] = cluster_records
    #
    # rep_names = []
    # for cluster_id in cluster_info.keys():
    #     contigname = ''
    #     for index, record in enumerate(cluster_info[cluster_id]):
    #         # representative record
    #         record = str(record)
    #         if record.endswith('... *') and record.__contains__('>Node_'):
    #             contigname = record.split('>')[1].replace('... *', '')
    #             contigname = contigname.split('-')[0]
    #     if len(cluster_info[cluster_id]) > 4 and contigname != '':
    #         rep_names.append(contigname)
    # print(rep_names)
    # print(len(rep_names))
    #
    # candidate_TE_path = tmp_output_dir + '/TE.consensus.fa'
    # with open(candidate_TE_path, 'w') as f_save:
    #     for name in rep_names:
    #         seq = new_candidate_contigs[name]
    #         f_save.write('>'+name+'\n'+seq+'\n')
    #
    # merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    # os.system('cat ' + candidate_TE_path + ' > ' + merge_pure)
    # ltr_retriever_seq = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna.mod.LTRlib.fa'
    # os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    #
    # tools_dir = os.getcwd() + '/tools'
    # cd_hit_command = tools_dir + '/cd-hit-est -s 0.8 -c 0.8 -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    # contigNames, contigs = read_fasta(merge_pure_consensus)
    # new_contigs = {k: v for k, v in sorted(contigs.items(), key=lambda item: len(item[1]))}
    # # store_fasta(new_contigs, home_dir+'/sort_lib/'+lib)
    # store_fasta(new_contigs, merge_pure_consensus)
    #
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    #
    # sample_name = 'oryza_sative'
    # reference = '/public/home/hpc194701009/repeat_detect_tools/WebTE/demo/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    # # sample_name = 'drerio'
    # # reference = '/public/home/hpc194701009/KmerRepFinder_test/genome/GCF_000002035.6_GRCz11_genomic.fna'
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + merge_pure_consensus + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(48) + ' -o ' + tmp_output_dir
    # os.system(TEClass_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step11. assign a family name for each classified TE consensus
    # classified_consensus_path = merge_pure_consensus + '.final.classified'
    #
    # # Step12. invoke RepeatMasker to align TE family to genome
    # RepeatMasker_Home = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker'
    # # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/drerio/drerio_test'
    # # RepeatMasker_output_dir = output_dir + '/drerio'
    # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/oryza_sative/oryza_sativa_test'
    # RepeatMasker_output_dir = output_dir + '/oryza_sative'
    # RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(48) \
    #                        + ' -lib ' + classified_consensus_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    # os.system('rm -rf ' + RepeatMasker_output_dir)
    # os.system(RepeatMasker_command)

    # sam_path_bwa = tmp_output_dir + '/repeats.filter_tandem.sam'
    # sam_paths = []
    # sam_paths.append(sam_path_bwa)
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds = get_alignment_info_v1(sam_paths, candidate_repeats_path)
    #
    # single_mapped_path = tmp_output_dir + '/repeats.merge.consensus.single.fa'
    # multiple_mapped_path = tmp_output_dir + '/repeats.merge.consensus.multiple.fa'
    # segmental_duplication_path = tmp_output_dir + '/segmental_duplication.fa'
    # merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # # single repeat probably be Chimeric repeat
    # with open(single_mapped_path, 'w') as f_save:
    #     for repeat_id in single_mapped_repeatIds:
    #         f_save.write('>' + repeat_id + '\n' + merge_repeat_contigs[repeat_id] + '\n')
    #
    # with open(multiple_mapped_path, 'w') as f_save:
    #     for repeat_id in multi_mapping_repeatIds:
    #         seq = merge_repeat_contigs[repeat_id]
    #         # seq = seq.replace('N', '')
    #         f_save.write('>' + repeat_id + '\n' + seq + '\n')
    #
    # with open(segmental_duplication_path, 'w') as f_save:
    #     for repeat_id in segmental_duplication_repeatIds:
    #         seq = merge_repeat_contigs[repeat_id]
    #         # seq = seq.replace('N', '')
    #         f_save.write('>' + repeat_id + '\n' + seq + '\n')
    # #
    # # # 06: merge
    # merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    # os.system('cat ' + multiple_mapped_path + ' > ' + merge_pure)
    # ltr_retriever_seq = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna.mod.LTRlib.fa'
    # os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    #
    # tools_dir = os.getcwd() + '/tools'
    # cd_hit_command = tools_dir + '/cd-hit-est -s 0.8 -c 0.8 -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    # contigNames, contigs = read_fasta(merge_pure_consensus)
    # new_contigs = {k: v for k, v in sorted(contigs.items(), key=lambda item: len(item[1]))}
    # # store_fasta(new_contigs, home_dir+'/sort_lib/'+lib)
    # store_fasta(new_contigs, merge_pure_consensus)
    #
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    #
    # sample_name = 'oryza_sative'
    # reference = '/public/home/hpc194701009/repeat_detect_tools/WebTE/demo/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    # # sample_name = 'drerio'
    # # reference = '/public/home/hpc194701009/KmerRepFinder_test/genome/GCF_000002035.6_GRCz11_genomic.fna'
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + merge_pure_consensus + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(48) + ' -o ' + tmp_output_dir
    # os.system(TEClass_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step11. assign a family name for each classified TE consensus
    # classified_consensus_path = merge_pure_consensus + '.final.classified'
    #
    # # Step12. invoke RepeatMasker to align TE family to genome
    # RepeatMasker_Home = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker'
    # # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/drerio/drerio_test'
    # # RepeatMasker_output_dir = output_dir + '/drerio'
    # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/oryza_sative/oryza_sativa_test'
    # RepeatMasker_output_dir = output_dir + '/oryza_sative'
    # RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(48) \
    #                        + ' -lib ' + classified_consensus_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    # os.system('rm -rf ' + RepeatMasker_output_dir)
    # os.system(RepeatMasker_command)