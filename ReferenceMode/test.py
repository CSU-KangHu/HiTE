import os

from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output, \
    get_multiple_alignment_repeat, split2cluster, cut_repeat_v1, judgeReduceThreads, get_ltr_suppl_from_ltrfinder, \
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v1

if __name__ == '__main__':
    curated_library = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/model_lib.fa'
    curated_contignames, curated_contigs = read_fasta(curated_library)
    output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output'
    # # extract from curated library
    model_lib_path = output_dir + '/model_lib.fa'
    with open(model_lib_path, 'w') as f_save:
        for name in curated_contignames:
            seq = curated_contigs[name]
            f_save.write('>' + name + '\n' + seq + '\n')

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