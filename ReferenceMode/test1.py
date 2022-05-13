import os

from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output, \
    get_multiple_alignment_repeat, split2cluster, cut_repeat_v1, judgeReduceThreads, get_ltr_suppl_from_ltrfinder, \
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v1

if __name__ == '__main__':
    # RepeatMasker_Home = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker'
    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/oryza_sative/CRD.2022-05-10.0-18-48'
    # repeats_path = tmp_output_dir + '/repeats.fa'
    # trf_dir = tmp_output_dir + '/trf_temp'
    # if not os.path.exists(trf_dir):
    #     os.makedirs(trf_dir)
    # (repeat_dir, repeat_filename) = os.path.split(repeats_path)
    # (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    # RepeatMasker_command = 'cd ' + trf_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(48) + ' -noint -x -html -gff -dir ' + trf_dir + ' ' + repeats_path
    # os.system(RepeatMasker_command)
    # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.masked'
    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(repeats_path)
    # repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # with open(repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('X')) / len(seq) < 0.5 and len(seq) >= 80:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')


    #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/drerio/CRD.2022-05-10.2-32-16'
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/oryza_sative/CRD.2022-05-11.11-45-41'
    candidate_repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa' 
    sam_path_bwa = tmp_output_dir + '/repeats.filter_tandem.sam' 
    sam_paths = []
    sam_paths.append(sam_path_bwa)
    unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds = get_alignment_info_v1(sam_paths, candidate_repeats_path)

    single_mapped_path = tmp_output_dir + '/repeats.merge.consensus.single.fa'
    multiple_mapped_path = tmp_output_dir + '/repeats.merge.consensus.multiple.fa'
    segmental_duplication_path = tmp_output_dir + '/segmental_duplication.fa'
    merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # single repeat probably be Chimeric repeat
    with open(single_mapped_path, 'w') as f_save:
        for repeat_id in single_mapped_repeatIds:
            f_save.write('>' + repeat_id + '\n' + merge_repeat_contigs[repeat_id] + '\n')

    with open(multiple_mapped_path, 'w') as f_save:
        for repeat_id in multi_mapping_repeatIds:
            seq = merge_repeat_contigs[repeat_id]
            # seq = seq.replace('N', '')
            f_save.write('>' + repeat_id + '\n' + seq + '\n')

    with open(segmental_duplication_path, 'w') as f_save:
        for repeat_id in segmental_duplication_repeatIds:
            seq = merge_repeat_contigs[repeat_id]
            # seq = seq.replace('N', '')
            f_save.write('>' + repeat_id + '\n' + seq + '\n')
    #
    # # 06: merge
    merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    os.system('cat ' + multiple_mapped_path + ' > ' + merge_pure)
    ltr_retriever_seq = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna.mod.LTRlib.fa'
    os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)

    tools_dir = os.getcwd() + '/tools'
    cd_hit_command = tools_dir + '/cd-hit-est -s 0.8 -c 0.8 -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)
    contigNames, contigs = read_fasta(merge_pure_consensus)
    new_contigs = {k: v for k, v in sorted(contigs.items(), key=lambda item: len(item[1]))}
    # store_fasta(new_contigs, home_dir+'/sort_lib/'+lib)
    store_fasta(new_contigs, merge_pure_consensus)

    merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'

    sample_name = 'oryza_sative'
    reference = '/public/home/hpc194701009/repeat_detect_tools/WebTE/demo/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    # sample_name = 'drerio'
    # reference = '/public/home/hpc194701009/KmerRepFinder_test/genome/GCF_000002035.6_GRCz11_genomic.fna'
    TEClass_home = os.getcwd() + '/classification'
    TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                      + ' --consensus ' + merge_pure_consensus + ' --genome ' + reference \
                      + ' --thread_num ' + str(48) + ' -o ' + tmp_output_dir
    os.system(TEClass_command)

    # --------------------------------------------------------------------------------------
    # Step11. assign a family name for each classified TE consensus
    classified_consensus_path = merge_pure_consensus + '.final.classified'

    # Step12. invoke RepeatMasker to align TE family to genome
    RepeatMasker_Home = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker'
    # output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/drerio/drerio_test'
    # RepeatMasker_output_dir = output_dir + '/drerio'
    output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/oryza_sative/oryza_sativa_test'
    RepeatMasker_output_dir = output_dir + '/oryza_sative'
    RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(48) \
                           + ' -lib ' + classified_consensus_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    os.system('rm -rf ' + RepeatMasker_output_dir)
    os.system(RepeatMasker_command)