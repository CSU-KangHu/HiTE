import argparse
import os
import time

from Util import read_fasta, Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output


def store_sequences(reference, candidate_pos, repeats_path,
                    max_complete_len, min_complete_len, direct2CompleteNames):
    refContigNames, refContigs = read_fasta(reference)

    # filter ltr/tir by min/max length
    candidate_info = {}
    for record in candidate_pos:
        query_name = record[0]
        target_name = record[1]
        start = int(record[2])
        end = int(record[5])
        seq_len = abs(end - start) + 1
        if seq_len > max_complete_len or seq_len < min_complete_len:
            continue
        if not candidate_info.__contains__(target_name):
            candidate_info[target_name] = []
        pos_array = candidate_info[target_name]
        pos_array.append((start, end, query_name))
        candidate_info[target_name] = pos_array

    # sort candidate_info by start and end position
    for target_name in candidate_info.keys():
        pos_array = candidate_info[target_name]
        pos_array.sort(key=lambda x: (x[0], x[1]))
        # merge redundant candidate sequence
        last = None
        unique_pos = {}
        for item in pos_array:
            if last is None:
                unique_pos[item] = 1
                last = item
            else:
                if item[0] >= last[0] and item[1] >= last[1]:
                    if item[0] == last[0]:
                        del unique_pos[last]
                    unique_pos[item] = 1
                    last = item
                elif item[1] < last[1]:
                    continue
        candidate_info[target_name] = list(unique_pos.keys())

    # step3. generate candidate complete repeats
    node_index = 0
    with open(repeats_path, 'w') as f_save:
        for target_name in candidate_info.keys():
            refContig = refContigs[target_name]
            for pos_info in candidate_info[target_name]:
                start = int(pos_info[0])
                end = int(pos_info[1])
                query_name = pos_info[2]
                sequence = refContig[start - 1: end]
                complete_name = 'Node_' + str(node_index) + '-len_' + str(len(sequence))
                if not sequence.__contains__('N'):
                    f_save.write('>' + complete_name + '\n' + sequence + '\n')
                    direct2CompleteNames[complete_name] = query_name
                    node_index += 1


def parse_blast_output(blastnResults_path, max_ltr_complete_len, min_ltr_complete_len,
                       max_tir_complete_len, min_tir_complete_len, min_non_ltr_length, max_non_ltr_length,
                       min_non_tir_length, max_non_tir_length, max_tir_direct_repeat_len,
                       min_identity, match_ratio,
                       reference, ltr_repeats_path, tir_repeats_path):

    # To facilite searching
    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}
    # step1. construct blast records clustering by query name
    query_records = {}
    with open(blastnResults_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            target_name = parts[1]
            identity = float(parts[2])
            match_base = int(parts[3])
            query_length = int(query_name.split('-')[1].split('_')[1])
            q_start = int(parts[6])
            q_end = int(parts[7])
            t_start = int(parts[8])
            t_end = int(parts[9])
            strand = '+' if t_end >= t_start else '-'
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            records = query_records[query_name]
            if not records.__contains__(target_name):
                records[target_name] = []
            same_target_records = records[target_name]
            if identity >= min_identity and float(match_base)/query_length >= match_ratio:
                item = (t_start, t_end, strand)
                if not same_target_records.__contains__(item):
                    same_target_records.append(item)
                    records[target_name] = same_target_records
                    query_records[query_name] = records

    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}
    # Node_0-len_5109 Node_0-len_5109 100.000 4651    0       0       459     5109    1       4651    0.0     8589
    # Node_0-len_5109 Node_30444-len_20481    100.000 217     0       0       1       217     20265   20481   1.37e-110       401

    # step2. find candidate ltr and tir sequence
    ltr_candidate_pos = []
    tir_candidate_pos = []
    for query_name in query_records.keys():
        query_length = int(query_name.split('-')[1].split('_')[1])
        records = query_records[query_name]
        for target_name in records.keys():
            same_chr_pos_records = records[target_name]
            for i in range(len(same_chr_pos_records)):
                for j in range(i + 1, len(same_chr_pos_records)):
                    # ltr direct repeats have the same direction
                    if same_chr_pos_records[j][2] == same_chr_pos_records[i][2]:
                        # filter cross
                        non_ltr_len = -1
                        if same_chr_pos_records[i][2] == '+':
                            if same_chr_pos_records[j][0] > same_chr_pos_records[i][1]:
                               non_ltr_len = same_chr_pos_records[j][0] - same_chr_pos_records[i][1]
                               if non_ltr_len >= min_non_ltr_length and non_ltr_len <= max_non_ltr_length:
                                   ltr_candidate_pos.append((query_name, target_name,
                                                             same_chr_pos_records[i][0], same_chr_pos_records[i][1],
                                                             same_chr_pos_records[j][0], same_chr_pos_records[j][1],
                                                             same_chr_pos_records[i][2]))
                            elif same_chr_pos_records[j][1] < same_chr_pos_records[i][0]:
                               non_ltr_len = same_chr_pos_records[i][0] - same_chr_pos_records[j][1]
                               if non_ltr_len >= min_non_ltr_length and non_ltr_len <= max_non_ltr_length:
                                   ltr_candidate_pos.append((query_name, target_name,
                                                             same_chr_pos_records[j][0], same_chr_pos_records[j][1],
                                                             same_chr_pos_records[i][0], same_chr_pos_records[i][1],
                                                             same_chr_pos_records[i][2]))
                        else:
                            if same_chr_pos_records[j][1] > same_chr_pos_records[i][0]:
                                non_ltr_len = same_chr_pos_records[j][1] - same_chr_pos_records[i][0]
                                if non_ltr_len >= min_non_ltr_length and non_ltr_len <= max_non_ltr_length:
                                    ltr_candidate_pos.append((query_name, target_name,
                                                              same_chr_pos_records[i][1], same_chr_pos_records[i][0],
                                                              same_chr_pos_records[j][1], same_chr_pos_records[j][0],
                                                              same_chr_pos_records[i][2]))
                            elif same_chr_pos_records[j][0] < same_chr_pos_records[i][1]:
                                non_ltr_len = same_chr_pos_records[i][1] - same_chr_pos_records[j][0]
                                if non_ltr_len >= min_non_ltr_length and non_ltr_len <= max_non_ltr_length:
                                    ltr_candidate_pos.append((query_name, target_name,
                                                              same_chr_pos_records[j][1], same_chr_pos_records[j][0],
                                                              same_chr_pos_records[i][1], same_chr_pos_records[i][0],
                                                              same_chr_pos_records[i][2]))
                    else:
                        # tir direct repeats have the reverse direction
                        # for tir flanking repeats commonly have several hundred bps
                        if query_length <= max_tir_direct_repeat_len:
                            if same_chr_pos_records[i][2] == '+':
                                if same_chr_pos_records[j][1] > same_chr_pos_records[i][1]:
                                    non_tir_len = same_chr_pos_records[j][1] - same_chr_pos_records[i][1]
                                    if non_tir_len >= min_non_tir_length and non_tir_len <= max_non_tir_length:
                                        tir_candidate_pos.append((query_name, target_name,
                                                                  same_chr_pos_records[i][0],
                                                                  same_chr_pos_records[i][1],
                                                                  same_chr_pos_records[j][1],
                                                                  same_chr_pos_records[j][0]))
                                elif same_chr_pos_records[j][0] < same_chr_pos_records[i][0]:
                                    non_tir_len = same_chr_pos_records[i][0] - same_chr_pos_records[j][0]
                                    if non_tir_len >= min_non_tir_length and non_tir_len <= max_non_tir_length:
                                        tir_candidate_pos.append((query_name, target_name,
                                                                  same_chr_pos_records[j][1],
                                                                  same_chr_pos_records[j][0],
                                                                  same_chr_pos_records[i][0],
                                                                  same_chr_pos_records[i][1]))
                            else:
                                if same_chr_pos_records[i][1] > same_chr_pos_records[j][1]:
                                    non_tir_len = same_chr_pos_records[i][1] - same_chr_pos_records[j][1]
                                    if non_tir_len >= min_non_tir_length and non_tir_len <= max_non_tir_length:
                                        tir_candidate_pos.append((query_name, target_name,
                                                                  same_chr_pos_records[j][0],
                                                                  same_chr_pos_records[j][1],
                                                                  same_chr_pos_records[i][1],
                                                                  same_chr_pos_records[i][0]))
                                elif same_chr_pos_records[i][0] < same_chr_pos_records[j][0]:
                                    non_tir_len = same_chr_pos_records[j][0] - same_chr_pos_records[i][0]
                                    if non_tir_len >= min_non_tir_length and non_tir_len <= max_non_tir_length:
                                        tir_candidate_pos.append((query_name, target_name,
                                                                  same_chr_pos_records[i][1],
                                                                  same_chr_pos_records[i][0],
                                                                  same_chr_pos_records[j][0],
                                                                  same_chr_pos_records[j][1]))
    # print(tir_candidate_pos)

    # ltr_direct_repeats_name 2 ltr_complete_repeats_name
    # it is obviously that ltr/tir_complete_repeats now include
    # False Positive sequence. To remove these sequences, we
    # need to record repeats names which construct these FP sequences.
    direct2CompleteNames = {}

    store_sequences(reference, ltr_candidate_pos, ltr_repeats_path,
                    max_ltr_complete_len, min_ltr_complete_len, direct2CompleteNames)

    store_sequences(reference, tir_candidate_pos, tir_repeats_path,
                    max_tir_complete_len, min_tir_complete_len, direct2CompleteNames)


    # candidate_info_path = '/public/home/hpc194701009/WebTE_Lib/Cicer/Cicer_arietinum_3827/kmerRepFinder/CRD.2022-01-19.15-32-32/ltr_tir.info'
    # with open(candidate_info_path, 'w') as f_save:
    #     for target_name in candidate_info.keys():
    #         for pos_info in candidate_info[target_name]:
    #             f_save.write(str(target_name)+'\t'+str(pos_info[0])+'\t'+str(pos_info[1])+'\t'+str(pos_info[2])+'\n')
    return direct2CompleteNames


def parse_blastx_output(blastxResults_path, protein_db_path,
                        domain_min_identity, domain_match_ratio):
    prot_contignames, prot_contigs = read_fasta(protein_db_path)
    # To facilite searching
    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}
    # step1. construct blast records clustering by query name
    query_records = {}
    with open(blastxResults_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            target_name = parts[1]
            target_name = target_name.replace(',', '')
            identity = float(parts[2])
            match_base = int(parts[3])
            target_length = len(prot_contigs[target_name])
            target_cover_ratio = float(match_base) / target_length
            q_start = int(parts[6])
            q_end = int(parts[7])
            t_start = int(parts[8])
            t_end = int(parts[9])
            strand = '+' if t_end >= t_start else '-'
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            records = query_records[query_name]
            if not records.__contains__(target_name):
                records[target_name] = []
            same_target_records = records[target_name]
            if identity >= domain_min_identity and target_cover_ratio >= domain_match_ratio:
                same_target_records.append((query_name, target_name, identity, target_cover_ratio,
                                            q_start, q_end, t_start, t_end, strand))
                records[target_name] = same_target_records
                query_records[query_name] = records

    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}

    # step2. find true positive complete LTR-RT/TIR sequence
    ltr_tir_domain_temp_map = {}
    ltr_tir_domain_map = {}
    for query_name in query_records.keys():
        records = query_records[query_name]
        if not ltr_tir_domain_temp_map.__contains__(query_name):
            ltr_tir_domain_temp_map[query_name] = []
        cur_ltr_domains = ltr_tir_domain_temp_map[query_name]
        for target_name in records.keys():
            # choose best match records, we consider both
            # identity and domain cover ratio
            same_target_records = records[target_name]
            best_match_scores = 0
            best_match_record = None
            for same_target_record in same_target_records:
                identity = same_target_record[2]
                target_cover_ratio = same_target_record[3]
                cur_score = identity * target_cover_ratio
                if cur_score > best_match_scores:
                    best_match_scores = cur_score
                    best_match_record = same_target_record
            if best_match_record is not None:
                cur_ltr_domains.append(best_match_record)
        # sort cur_ltr_domains by query start
        cur_ltr_domains.sort(key=lambda x: (x[4], x[5]))
        if len(cur_ltr_domains) > 0:
            ltr_tir_domain_map[query_name] = cur_ltr_domains

    return ltr_tir_domain_map


if __name__ == '__main__':
    log = Logger('kmerRepFinder.log', level='debug')

    # 1.parse args
    parser = argparse.ArgumentParser(description='run TEFinder...')
    parser.add_argument('--min_TE_len', metavar='Minimum length of TE sequences',
                        help='input minimum length of TE sequences')
    parser.add_argument('--min_ltr_complete_len', metavar='Minimum length of the complete LTR-RT sequences',
                        help='input minimum length of the complete LTR-RT sequences')
    parser.add_argument('--max_ltr_complete_len', metavar='Maximum  length of the complete LTR-RT sequences',
                        help='input maximum  length of the complete LTR-RT sequences')
    parser.add_argument('--min_ltr_direct_repeat_len', metavar='Minimum length of LTR direct repeat',
                        help='input minimum length of LTR direct repeat')
    parser.add_argument('--max_ltr_direct_repeat_len', metavar='Maximum length of LTR direct repeat',
                        help='input maximum length of LTR direct repeat')

    parser.add_argument('--min_tir_complete_len', metavar='Minimum length of the complete TIR sequences',
                        help='input minimum length of the complete TIR sequences')
    parser.add_argument('--max_tir_complete_len', metavar='Maximum  length of the complete TIR sequences',
                        help='input maximum  length of the complete TIR sequences')
    parser.add_argument('--min_tir_direct_repeat_len', metavar='Minimum length of TIR direct repeat',
                        help='input minimum length of TIR direct repeat')
    parser.add_argument('--max_tir_direct_repeat_len', metavar='Maximum length of TIR direct repeat',
                        help='input maximum length of TIR direct repeat')

    parser.add_argument('--tmp_output_dir', metavar='Temp output directory',
                        help='input temp output directory')
    parser.add_argument('-R', metavar='Reference Path',
                        help='input reference path')
    parser.add_argument('-t', metavar='Thread num',
                        help='input thread num')

    parser.add_argument('--domain_min_identity', metavar='Minimum identity of True Positive LTR and TIR',
                        help='input minimum identity of True Positive LTR and TIR')
    parser.add_argument('--domain_match_ratio', metavar='Minimum cover ratio of domain for True Positive LTR and TIR',
                        help='input minimum cover ratio of domain for True Positive LTR and TIR')
    parser.add_argument('-s', metavar='sensitive mode',
                        help='sensitive mode, default 0')
    parser.add_argument('--long_repeat_threshold', metavar='long_repeat_threshold',
                        help='long_repeat_threshold, default 2000')
    parser.add_argument('--tools_dir', metavar='tools dir',
                        help='input tools dir')



    args = parser.parse_args()

    min_TE_len = int(args.min_TE_len)

    # params get from LTRDetector
    min_ltr_complete_len = int(args.min_ltr_complete_len)
    max_ltr_complete_len = int(args.max_ltr_complete_len)
    min_ltr_direct_repeat_len = int(args.min_ltr_direct_repeat_len)
    max_ltr_direct_repeat_len = int(args.max_ltr_direct_repeat_len)

    # Class II transposons range in length from 1,000 to as many as 40,000 base pairs.
    # https://www.britannica.com/science/transposon
    min_tir_complete_len = int(args.min_tir_complete_len)
    max_tir_complete_len = int(args.max_tir_complete_len)
    min_tir_direct_repeat_len = int(args.min_tir_direct_repeat_len)
    max_tir_direct_repeat_len = int(args.max_tir_direct_repeat_len)

    tmp_output_dir = args.tmp_output_dir
    long_repeat_threshold = int(args.long_repeat_threshold)
    tools_dir = args.tools_dir
    sensitive_mode = args.s

    is_sensitive = False
    if sensitive_mode == '1':
        is_sensitive = True
    # --------------------------------------------------------------------------------------
    # Step1. get LTR direct repeats from repeats consensus
    # len >= 100 and len <= 6000
    consensus_path = tmp_output_dir + '/repeats_consensus.fa'
    direct_candidate = tmp_output_dir + '/direct_repeat.candidate.fa'
    consensusNames, consensusContigs = read_fasta(consensus_path)
    with open(direct_candidate, 'w') as f_save:
        for consensusName in consensusNames:
            seq = consensusContigs[consensusName]
            if len(seq) >= min_ltr_direct_repeat_len and len(seq) <= max_ltr_direct_repeat_len:
                f_save.write('>'+consensusName+'\n'+seq+'\n')

    reference = args.R

    ltr_repeats_path = tmp_output_dir + '/ltr_repeats.fa'
    tir_repeats_path = tmp_output_dir + '/tir_repeats.fa'
    ltr_tir_repeats_consensus_path = tmp_output_dir + '/ltr_tir_repeats.consensus.fa'
    blast_program_dir = os.getcwd() + '/tools/rmblast-2.9.0-p2'
    blastnResults_path = tmp_output_dir + '/tmpBlastResults.ltr_tir.out'
    threads = args.t

    # --------------------------------------------------------------------------------------
    # Step2. align ltr direct repeats to reference
    # get candidate ltr complete sequence (including protein domain)
    starttime = time.time()
    if is_sensitive:
        log.logger.debug('Start TEFinder module1: use blast to execute repeats consensus align to reference')
        makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + reference
        align_command = blast_program_dir + '/bin/blastn -db ' + reference + ' -num_threads ' + \
                        threads + ' -query ' + direct_candidate + ' -outfmt 6 > ' + blastnResults_path
        log.logger.debug(makedb_command)
        os.system(makedb_command)
        log.logger.debug(align_command)
        os.system(align_command)
    else:
        log.logger.debug('TEFinder module1: use minimap2+bwa to execute repeats consensus align to reference')
        direct_candidate_minimap2 = tmp_output_dir + '/direct_repeat.candidate.minimap2.fa'
        direct_candidate_bwa = tmp_output_dir + '/direct_repeat.candidate.bwa.fa'
        split_repeats(direct_candidate, long_repeat_threshold, direct_candidate_minimap2, direct_candidate_bwa)
        use_align_tools = 'minimap2'
        sam_path_minimap2 = run_alignment(direct_candidate_minimap2, reference, use_align_tools, threads, tools_dir)
        use_align_tools = 'bwa'
        sam_path_bwa = run_alignment(direct_candidate_bwa, reference, use_align_tools, threads, tools_dir)
        sam_paths = []
        sam_paths.append((sam_path_minimap2, direct_candidate_minimap2))
        sam_paths.append((sam_path_bwa, direct_candidate_bwa))
        not_multi_mapping_repeatIds = []
        generate_blastlike_output(sam_paths, blastnResults_path, not_multi_mapping_repeatIds)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("TEFinder module1: Blastn alignment running time: %.8s s" % (dtime))

    # parse blastn output, find LTR complete sequence
    min_non_ltr_length = 0
    max_non_ltr_length = max_ltr_complete_len - 2 * min_ltr_direct_repeat_len

    min_non_tir_length = 0
    max_non_tir_length = max_tir_complete_len - 2 * min_tir_direct_repeat_len

    # min_non_ltr_length = min_ltr_complete_len - 2*min_ltr_direct_repeat_len
    # max_non_ltr_length = max_ltr_complete_len - 2*max_ltr_direct_repeat_len
    #
    # min_non_tir_length = min_tir_complete_len - 2 * min_tir_direct_repeat_len
    # max_non_tir_length = max_tir_complete_len - 2 * max_tir_direct_repeat_len

    # I think repeats consensus should match high
    # enough to ensure this repeats is confident. at least 95% confidence
    # min_identity = 95
    # match_ratio = 0.95
    min_identity = 90
    match_ratio = 0.90
    # direct2CompleteNames is the map information about repeats name to complete ltr name
    direct2CompleteNames = parse_blast_output(blastnResults_path, max_ltr_complete_len, min_ltr_complete_len,
                                              max_tir_complete_len, min_tir_complete_len, min_non_ltr_length, max_non_ltr_length,
                                              min_non_tir_length, max_non_tir_length, max_tir_direct_repeat_len,
                                              min_identity, match_ratio,
                                              reference, ltr_repeats_path, tir_repeats_path)

    # --------------------------------------------------------------------------------------
    # Step3. get candidate TE sequences
    # merge_repeat_sequences = tmp_output_dir + '/repeats.merge.fa'
    # merge_command1 = 'cat ' + ltr_repeats_path + ' >> ' + merge_repeat_sequences
    # merge_command2 = 'cat ' + tir_repeats_path + ' >> ' + merge_repeat_sequences
    # merge_command3 = 'cat ' + consensus_path + ' >> ' + merge_repeat_sequences
    # os.system(merge_command1)
    # os.system(merge_command2)
    # os.system(merge_command3)
