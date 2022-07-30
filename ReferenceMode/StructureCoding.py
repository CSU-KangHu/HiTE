import os

from Util import read_fasta, store_fasta

# LTR range from 100-1000
# complete LTR-RTs range from 1100-16000

def print_seqs(header, sequence, length, outfile):
    print('>' + header, file=outfile)
    while len(sequence) > 0:
        print(sequence[:length], file=outfile)
        sequence = sequence[length:]

if __name__ == '__main__':
    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-07-10.21-5-55'
    tmp_output_dir = '/home/hukang/KRF_output/dmel/CRD.2022-07-09.16-8-30'
    TE_merge_path = tmp_output_dir + '/TE.merge.fa'
    TE_merge_multi_line_path = tmp_output_dir + '/TE.merge.ml.fa'
    contigNames, contigs = read_fasta(TE_merge_path)
    outfile = open(TE_merge_multi_line_path, 'w')  # open outfile for writing
    for name in contigNames:
        print_seqs(name, contigs[name], 50, outfile)


    # 1. run TRsearch
    TRsearch_dir = '/home/hukang/repeat_detect_tools/REPET_linux-x64-3.0/bin'
    TRsearch_command = TRsearch_dir + '/TRsearch ' + TE_merge_multi_line_path
    #os.system(TRsearch_command)


    # 2.parse TRsearch results
    TR_out = tmp_output_dir + '/TE.merge.ml.fa.TR.set'

    reference = '/home/hukang/KmerRepFinder_git/KmerRepFinder/demo/dmel-all-chromosome-r5.43.fasta'

    # ref_contigNames, ref_contigs = read_fasta(reference)
    #-----------------------------get TE sequences which have domain--------------------------------------
    TE_merge_classified_path = TE_merge_multi_line_path + '.final.classified'
    sample_name = 'dmel'
    TEClass_home = os.getcwd() + '/classification'
    TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                      + ' --consensus ' + TE_merge_multi_line_path + ' --genome ' + reference \
                      + ' --thread_num ' + str(40) + ' -o ' + tmp_output_dir
    #os.system(TEClass_command)

    query_sets = set()
    final_tmpBlastX_path = TE_merge_multi_line_path + '.tmpBlastXResults.out.bxsummary'

    with open(final_tmpBlastX_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            identity_parts = parts[6].split('/')
            identity = float(identity_parts[0])/ int(identity_parts[1])
            if identity >= 0.8:
                query_sets.add(query_name)

    ltr_count = 0
    TE_merge_contigNames, TE_merge_contigs = read_fasta(TE_merge_multi_line_path)
    TE_domain_path = tmp_output_dir + '/TE.domain.fa'
    TE_domain_contigs = {}
    for query_name in query_sets:
        TE_domain_contigs[query_name] = TE_merge_contigs[query_name]
        if query_name.__contains__('ltr'):
            ltr_count +=1
    store_fasta(TE_domain_contigs, TE_domain_path)
    print('ltr_count: %d' %ltr_count)

    # # -----------------------------find sequences with obvious LTR/TIR structure--------------------------------------
    # TE_merge_contigNames, TE_merge_contigs = read_fasta(TE_merge_multi_line_path)
    #
    # type_sets = set()
    # # query_structs = {query_name: {repeat_id: [r1, r2]} }
    # LTR_query_structs = {}
    # with open(TR_out, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         repeat_id = parts[0]
    #         query_name = parts[2]
    #         query_len = int(query_name.split('-len_')[1])
    #         start_pos = int(parts[3])
    #         end_pos = int(parts[4])
    #         struct_info = parts[1]
    #         struct_parts = struct_info.split('|')
    #         s_type = struct_parts[0].replace('Rep'+str(repeat_id), '')
    #         s_len = int(struct_parts[1].replace('len=', ''))
    #         s_identity = float(struct_parts[2].replace('id=', ''))
    #         TE_len = int(struct_parts[3].replace('lenTE=', ''))
    #         struct_tuple = (s_type, s_len, s_identity, TE_len, start_pos, end_pos)
    #         type_sets.add(s_type)
    #
    #         if float(struct_tuple[3])/query_len >= 0.95:
    #             if (struct_tuple[0] == 'LTR' or struct_tuple[0] == 'termLTR') \
    #                     and (struct_tuple[1] >= 100 and struct_tuple[1] <= 5000):
    #                 if not LTR_query_structs.__contains__(query_name):
    #                     LTR_query_structs[query_name] = {}
    #                 struct_records = LTR_query_structs[query_name]
    #                 if not struct_records.__contains__(repeat_id):
    #                     struct_records[repeat_id] = []
    #                 struct_tuples = struct_records[repeat_id]
    #                 struct_tuples.append(struct_tuple)
    #
    # #print(type_sets)
    # # choose the most probably LTR (longer in LTR len)
    # best_LTR_query_structs = {}
    # for query_name in LTR_query_structs.keys():
    #     struct_records = LTR_query_structs[query_name]
    #     query_len = int(query_name.split('-len_')[1])
    #     if len(struct_records) <= 1:
    #         best_LTR_query_structs[query_name] = list(struct_records.values())[0]
    #     else:
    #         best_record = None
    #         for repeat_id in struct_records.keys():
    #             struct_tuples = struct_records[repeat_id]
    #             first_tuple = struct_tuples[0]
    #             second_tuple = struct_tuples[1]
    #             if best_record is None:
    #                 best_record = struct_tuples
    #             elif first_tuple[1] > best_record[0][1]:
    #                 best_record = struct_tuples
    #         best_LTR_query_structs[query_name] = best_record
    #
    # TE_ltr_path = tmp_output_dir + '/TE.domain_ltr.fa'
    # for query_name in best_LTR_query_structs:
    #     TE_domain_contigs[query_name] = TE_merge_contigs[query_name]
    # store_fasta(TE_domain_contigs, TE_ltr_path)




    # #-----------------------------TR structure analyze--------------------------------------
    # TE_merge_contigNames, TE_merge_contigs = read_fasta(TE_merge_path)
    #
    # type_sets = set()
    # # query_structs = {query_name: {repeat_id: [r1, r2]} }
    # LTR_query_structs = {}
    # TIR_query_structs = {}
    # with open(TR_out, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         repeat_id = parts[0]
    #         query_name = parts[2]
    #         query_len = int(query_name.split('-len_')[1])
    #         start_pos = int(parts[3])
    #         end_pos = int(parts[4])
    #         struct_info = parts[1]
    #         struct_parts = struct_info.split('|')
    #         s_type = struct_parts[0].replace('Rep'+str(repeat_id), '')
    #         s_len = int(struct_parts[1].replace('len=', ''))
    #         s_identity = float(struct_parts[2].replace('id=', ''))
    #         TE_len = int(struct_parts[3].replace('lenTE=', ''))
    #         struct_tuple = (s_type, s_len, s_identity, TE_len, start_pos, end_pos)
    #         type_sets.add(s_type)
    #
    #         if float(struct_tuple[3])/query_len >= 0.95:
    #             if (struct_tuple[0] == 'LTR' or struct_tuple[0] == 'termLTR') \
    #                     and (struct_tuple[1] >= 100 and struct_tuple[1] <= 1000):
    #                     #and (query_len >= 1100 and query_len <= 16000):
    #                 if not LTR_query_structs.__contains__(query_name):
    #                     LTR_query_structs[query_name] = {}
    #                 struct_records = LTR_query_structs[query_name]
    #                 if not struct_records.__contains__(repeat_id):
    #                     struct_records[repeat_id] = []
    #                 struct_tuples = struct_records[repeat_id]
    #                 struct_tuples.append(struct_tuple)
    #             elif (struct_tuple[0] == 'termTIR' or struct_tuple[0] == 'non-termTIR'):
    #                 if not TIR_query_structs.__contains__(query_name):
    #                     TIR_query_structs[query_name] = {}
    #                 struct_records = TIR_query_structs[query_name]
    #                 if not struct_records.__contains__(repeat_id):
    #                     struct_records[repeat_id] = []
    #                 struct_tuples = struct_records[repeat_id]
    #                 struct_tuples.append(struct_tuple)
    #
    # #print(type_sets)
    # # choose the most probably LTR (longer in LTR len)
    # best_LTR_query_structs = {}
    # for query_name in LTR_query_structs.keys():
    #     struct_records = LTR_query_structs[query_name]
    #     query_len = int(query_name.split('-len_')[1])
    #     if len(struct_records) <= 1:
    #         best_LTR_query_structs[query_name] = struct_records.values()[0]
    #     else:
    #         best_record = None
    #         for repeat_id in struct_records.keys():
    #             struct_tuples = struct_records[repeat_id]
    #             first_tuple = struct_tuples[0]
    #             second_tuple = struct_tuples[1]
    #             if best_record is None:
    #                 best_record = struct_tuples
    #             elif first_tuple[1] > best_record[0][1]:
    #                 best_record = struct_tuples
    #         best_LTR_query_structs[query_name] = best_record
    #
    # # choose the most probably TIR (longer in TE len)
    # best_TIR_query_structs = {}
    # for query_name in TIR_query_structs.keys():
    #     if best_LTR_query_structs.__contains__(query_name):
    #         continue
    #     struct_records = TIR_query_structs[query_name]
    #     query_len = int(query_name.split('-len_')[1])
    #     if len(struct_records) <= 1:
    #         best_TIR_query_structs[query_name] = struct_records.values()[0]
    #     else:
    #         best_record = None
    #         for repeat_id in struct_records.keys():
    #             struct_tuples = struct_records[repeat_id]
    #             first_tuple = struct_tuples[0]
    #             second_tuple = struct_tuples[1]
    #             if best_record is None:
    #                 best_record = struct_tuples
    #             elif first_tuple[3] > best_record[0][3]:
    #                 best_record = struct_tuples
    #         best_TIR_query_structs[query_name] = best_record
    # print(len(best_TIR_query_structs))
    #
    # #node_index = 0
    # LTR_no_struct = {}
    # clear_LTR_seqs = {}
    # for query_name in best_LTR_query_structs:
    #     best_record = best_LTR_query_structs[query_name]
    #     first_tuple = best_record[0]
    #     second_tuple = best_record[1]
    #
    #     leftLTR_start = first_tuple[4]
    #     leftLTR_end = first_tuple[5]
    #
    #     ltr_len = int(first_tuple[1])
    #     ltr_identity = float(first_tuple[2])
    #     TE_len = int(first_tuple[3])
    #
    #     rightLTR_start = second_tuple[4]
    #     rightLTR_end = second_tuple[5]
    #     query_seq = TE_merge_contigs[query_name]
    #     clear_LTR_seq = query_seq[leftLTR_start-1: rightLTR_end]
    #     #print(, TE_len)
    #     # LTR_info = (left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end, complete_LTR-RT_len, dna_sequence)
    #     LTR_info = (1, ltr_len, len(clear_LTR_seq)-ltr_len+1, len(clear_LTR_seq), len(clear_LTR_seq), clear_LTR_seq)
    #     # clear_LTR_seqs['N_'+str(node_index)+'-len_'+str(len(clear_LTR_seq))] = clear_LTR_seq
    #     clear_LTR_seqs[query_name] = clear_LTR_seq
    #     LTR_no_struct[query_name] = query_seq[leftLTR_end: rightLTR_start-1]
    #
    # TIR_no_struct = {}
    # clear_TIR_seqs = {}
    # for query_name in best_TIR_query_structs:
    #     best_record = best_TIR_query_structs[query_name]
    #     first_tuple = best_record[0]
    #     second_tuple = best_record[1]
    #
    #     leftTIR_start = first_tuple[4]
    #     leftTIR_end = first_tuple[5]
    #
    #     tir_len = int(first_tuple[1])
    #     tir_identity = float(first_tuple[2])
    #     TE_len = int(first_tuple[3])
    #
    #     rightTIR_start = second_tuple[4]
    #     rightTIR_end = second_tuple[5]
    #     query_seq = TE_merge_contigs[query_name]
    #     clear_TIR_seq = query_seq[leftTIR_end - 1: rightTIR_start]
    #     # print(, TE_len)
    #     # TIR_info = (left_tir_start, left_tir_end, right_tir_start, right_tir_end, complete_TIR_len, dna_sequence)
    #     TIR_info = (1, tir_len, len(clear_TIR_seq) - tir_len + 1, len(clear_TIR_seq), len(clear_TIR_seq), clear_TIR_seq)
    #     # clear_LTR_seqs['N_'+str(node_index)+'-len_'+str(len(clear_LTR_seq))] = clear_LTR_seq
    #     clear_TIR_seqs[query_name] = clear_TIR_seq
    #     TIR_no_struct[query_name] = query_seq[leftTIR_start: rightTIR_end - 1]
    #
    # LTR_no_struct_path = tmp_output_dir + '/LTR.no_struct.fa'
    # store_fasta(LTR_no_struct, LTR_no_struct_path)
    #
    # TIR_no_struct_path = tmp_output_dir + '/TIR.no_struct.fa'
    # store_fasta(TIR_no_struct, TIR_no_struct_path)
    #
    # # --------------------------------------------------------------------------------------
    # # Step7. run TE classification to classify TE family
    # sample_name = 'dmel'
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + LTR_no_struct_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(40) + ' -o ' + tmp_output_dir
    # os.system(TEClass_command)
    #
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + TIR_no_struct_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(40) + ' -o ' + tmp_output_dir
    # os.system(TEClass_command)
    #
    # LTR_no_struct_classified_path = LTR_no_struct_path + '.final.classified'
    # TIR_no_struct_classified_path = TIR_no_struct_path + '.final.classified'
    # ltr_classified_contigNames, ltr_classified_contigs = read_fasta(LTR_no_struct_classified_path)
    # tir_classified_contigNames, tir_classified_contigs = read_fasta(TIR_no_struct_classified_path)
    #
    # TE_complete_contigs = {}
    # for name in ltr_classified_contigNames:
    #     if name.__contains__('LTR'):
    #         query_name = name.split('#')[0]
    #         TE_complete_contigs[query_name] = clear_LTR_seqs[query_name]
    # for name in tir_classified_contigNames:
    #     if name.__contains__('DNA'):
    #         query_name = name.split('#')[0]
    #         TE_complete_contigs[query_name] = clear_TIR_seqs[query_name]
    #
    #
    # complete_TE_path = tmp_output_dir + '/TE.complete.fa'
    # store_fasta(TE_complete_contigs, complete_TE_path)
    #
    #
    # TE_merge_classified_path = TE_merge_path + '.final.classified'
    # sample_name = 'dmel'
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + TE_merge_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(40) + ' -o ' + tmp_output_dir
    # os.system(TEClass_command)
    #
    # TE_classified_path = tmp_output_dir + '/TE.classified.fa'
    # TE_merge_classified_contigNames, TE_merge_classified_contigs = read_fasta(TE_merge_classified_path)
    #
    # with open(TE_classified_path, 'w') as f_save:
    #     for name in TE_merge_classified_contigNames:
    #         parts = name.split('#')
    #         query_name = parts[0]
    #         class_name = parts[1]
    #         if TE_complete_contigs.__contains__(query_name):
    #             sequence = TE_complete_contigs[query_name]
    #             f_save.write('>' + name + '\n' + sequence + '\n')
    #         elif class_name != 'Unknown':
    #             sequence = TE_merge_classified_contigs[name]
    #             f_save.write('>' + name + '\n' + sequence + '\n')


