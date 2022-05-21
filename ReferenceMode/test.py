import codecs
import json
import os
import sys
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


    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-20.19-53-55'
    reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-20.19-53-55'
    reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'
    output_dir = tmp_output_dir
    repeats_consensus = tmp_output_dir + '/repeats.consensus.fa'
    skip_threshold = 100
    identity_threshold = 0.90
    length_similarity_cutoff = 0.90
    partitions_num = 48
    # try to chain all fragments
    # --------------------------------------------------------------------------------------
    # newly strategy: 2022-04-29 by Kang Hu
    # 01: use bwa to get single mapped sequence
    candidate_repeats_path = repeats_consensus
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    use_align_tools = 'bwa'
    tools_dir = os.getcwd() + '/tools'
    sam_path_bwa = run_alignment(candidate_repeats_path, reference, use_align_tools, 48, tools_dir)
    sam_paths = []
    sam_paths.append(sam_path_bwa)
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info(sam_paths)
    new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, candidate_repeats_path)
    repeat_freq_path = tmp_output_dir + '/repeats.freq.fa'
    merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # single repeat probably be Chimeric repeat
    with open(repeat_freq_path, 'w') as f_save:
        for repeat_id in new_mapping_repeatIds.keys():
            freq = new_mapping_repeatIds[repeat_id][0]
            seq = merge_repeat_contigs[repeat_id]
            if freq <= 1 or (freq < 5 and len(seq) < 80):
                continue
            f_save.write('>' + repeat_id + '\tcopies=' + str(freq + 1) + '\n' + seq + '\n')

    # --------------------------------------------------------------------------------------
    # skip variation between fragments: 2022-05-19 by Kang Hu
    # sort by start and end position
    sorted_fragments = {}
    for ref_name in query_position.keys():
        position_array = query_position[ref_name]
        position_array.sort(key=lambda x: (x[1], x[2]))
        sorted_fragments[ref_name] = position_array

    # find all regions could be connected, threshold = 200bp
    # region_list keeps all regions, which include all fragments can be connected
    # region_list = {
    # R1: {
    # ref_name: [(F1, start, end),(F2, start, end),(F3, start, end),(F4, start, end)]
    # }
    # }
    region_list = {}
    last_end_pos = -1
    region_index = 0
    for ref_name in sorted_fragments.keys():
        cur_ref_fragments = sorted_fragments[ref_name]
        for item in cur_ref_fragments:
            query_name = item[0]
            start = item[1]
            end = item[2]
            region_id = 'R' + str(region_index)
            if not region_list.__contains__(region_id):
                region_list[region_id] = {}
            cur_region_dict = region_list[region_id]
            if not cur_region_dict.__contains__(ref_name):
                cur_region_dict[ref_name] = []
            cur_region_list = cur_region_dict[ref_name]
            if last_end_pos == -1:
                # first fragment add into cur_region_list directly
                cur_region_list.append((query_name, start, end))
            else:
                # cur fragment close to last fragment
                if start - last_end_pos < skip_threshold:
                    cur_region_list.append((query_name, start, end))
                else:
                    # cur fragment far from last fragment, start a new region
                    region_index += 1
                    region_id = 'R' + str(region_index)
                    if not region_list.__contains__(region_id):
                        region_list[region_id] = {}
                    cur_region_dict = region_list[region_id]
                    if not cur_region_dict.__contains__(ref_name):
                        cur_region_dict[ref_name] = []
                    cur_region_list = cur_region_dict[ref_name]
                    cur_region_list.append((query_name, start, end))
                    cur_region_dict[ref_name] = cur_region_list
                    region_list[region_id] = cur_region_dict
            cur_region_dict[ref_name] = cur_region_list
            region_list[region_id] = cur_region_dict
            last_end_pos = end
    print('finish region_list...')
    print(len(region_list))

    # start parallelization to generate fragments combination in each region
    # region_cluster = split2cluster(list(region_list.items()), partitions_num)
    # # region_combination keeps all combination of fragment in one region
    # # e.g., region_combination = {
    # # R1: {
    # # max_combination_len : 3,
    # # combinations: {
    # # c=3: [F1F2F3],
    # # c=2: [F1F2, F2F3],
    # # c=1: [F1,F2,F3]
    # # }
    # # }
    # # }
    #
    # # frag_hash keeps all fragment combination information: combination_len, reference, start, end
    # # e.g., frag_hash = {
    # # F1F2: {
    # # R1: (c=2, ref_name, start, end)
    # # R2: (c=2, ref_name, start, end)
    # # }
    # # }
    region_combination_tmp = tmp_output_dir + '/region_combination_tmp'
    if not os.path.exists(region_combination_tmp):
        os.makedirs(region_combination_tmp)
    MAX_JOBS_IN_QUEUE = 500
    #split_region_combination_size = 1024 * 1024 * 1024 # (1G)
    split_region_combination_size = 1024  # (1K)
    region_combination = {}
    region_combination_index = 0
    frag_hash = {}
    ex = ProcessPoolExecutor(partitions_num)
    total_region_list = list(region_list.items())
    jobs_left = len(total_region_list)
    jobs_iter = iter(total_region_list)
    jobs = {}
    while jobs_left:
        for region_item in jobs_iter:
            job = ex.submit(getRegionCombination, region_item)
            jobs[job] = 1
            if len(jobs) > MAX_JOBS_IN_QUEUE:
                break  # limit the job submission for now job

        for job in as_completed(jobs):
            jobs_left -= 1
            part_region_combination, part_frag_hash = job.result()
            region_combination.update(part_region_combination)
            for combine_name in part_frag_hash.keys():
                part_dict = part_frag_hash[combine_name]
                if not frag_hash.__contains__(combine_name):
                    frag_hash[combine_name] = part_dict
                else:
                    exist_dict = frag_hash[combine_name]
                    exist_dict.update(part_dict)
                    frag_hash[combine_name] = exist_dict
            del jobs[job]

            # store region_combination into file when over 1G
            region_combination_size = float(sys.getsizeof(region_combination))
            if region_combination_size >= split_region_combination_size:
                region_combination_file = region_combination_tmp + '/rc_' + str(region_combination_index) + '.csv'
                with codecs.open(region_combination_file, 'w', encoding='utf-8') as f:
                    json.dump(region_combination, f)
                region_combination_index += 1
                region_combination.clear()
    if len(region_combination) > 0:
        region_combination_file = region_combination_tmp + '/rc_' + str(region_combination_index) + '.csv'
        with codecs.open(region_combination_file, 'w', encoding='utf-8') as f:
            json.dump(region_combination, f)
        region_combination_index += 1
        region_combination.clear()
    ex.shutdown(wait=True)

    # start parallelization to get candidate combine fragments
    regionContigs = {}
    refNames, refContigs = read_fasta(reference)
    ex = ProcessPoolExecutor(partitions_num)
    jobs = {}
    partiton_index = 0
    for rc_index in range(region_combination_index):
        region_combination_file = region_combination_tmp + '/rc_' + str(rc_index) + '.csv'
        file = open(region_combination_file, 'r')
        js = file.read()
        region_combination = json.loads(js)

        # region_combination_cluster = split2cluster(list(region_combination.items()), partitions_num)
        total_region_combination = list(region_combination.items())
        jobs_left = len(total_region_combination)
        jobs_iter = iter(total_region_combination)

        while jobs_left:
            for region_combination_item in jobs_iter:
                # create part variable to reduce memory copy
                part_frag_hash = {}
                region_id = region_combination_item[0]
                cur_region_combination = region_combination_item[1]
                combinations = cur_region_combination['combinations']
                max_combination_len = cur_region_combination['max_combination_len']
                for c in range(max_combination_len, 0, -1):
                    cur_combinations = combinations[str(c)]
                    for combine_name in cur_combinations:
                        part_frag_hash[combine_name] = frag_hash[combine_name]

                region_dict = region_list[region_id]
                # ref_name = list(region_dict.keys())[0]
                # ref_seq = refContigs[ref_name]
                # submit job
                job = ex.submit(getCombineFragments, region_combination_item, part_frag_hash, identity_threshold,
                                length_similarity_cutoff, refContigs, region_dict, tmp_output_dir, blast_program_dir,
                                partiton_index)
                jobs[job] = 1
                partiton_index += 1
                if len(jobs) > MAX_JOBS_IN_QUEUE:
                    break  # limit the job submission for now job

            for job in as_completed(jobs):
                jobs_left -= 1
                part_regionContigs = job.result()
                del jobs[job]
                regionContigs.update(part_regionContigs)

    ex.shutdown(wait=True)

    connected_repeats = tmp_output_dir + '/repeats.connect.fa'
    with open(connected_repeats, 'w') as f_save:
        node_index = 0
        for region_id in regionContigs.keys():
            for frag in regionContigs[region_id]:
                seq = refContigs[frag[1]][frag[2]: frag[3]]
                f_save.write('>Node_' + str(node_index) + '\n' + seq + '\n')
                node_index += 1

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