import argparse
import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

def generate_full_length_out_v2(alignment, TE_lib, reference, tmp_output_dir, full_length_threshold, category, type_out="RM"):
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    threads = 1
    full_length_annotations, copies_direct, all_query_copies = get_full_length_copies_from_blastn_v2(TE_lib, reference,
                                                                                                     alignment,
                                                                                                     tmp_output_dir,
                                                                                                     threads, type_out,
                                                                                                     full_length_threshold,
                                                                                                     category)
    lines = set()
    for query_name in full_length_annotations.keys():
        query_name = str(query_name)
        for copy_annotation in full_length_annotations[query_name]:
            chr_pos = copy_annotation[0]
            annotation = copy_annotation[1]
            parts = chr_pos.split(':')
            chr_name = parts[0]
            chr_pos_parts = parts[1].split('-')
            chr_start = int(chr_pos_parts[0]) + 1
            chr_end = int(chr_pos_parts[1])
            new_line = (query_name, chr_name, chr_start, chr_end, -1)
            lines.add(new_line)
    return lines

def parse_blastn_out(blastOut):
    query_records = {}
    with open(blastOut, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            info_parts = line.split('\t')
            query_name = info_parts[0].split('#')[0]
            subject_name = info_parts[1]
            q_start = int(info_parts[6])
            q_end = int(info_parts[7])
            s_start = int(info_parts[8])
            s_end = int(info_parts[9])
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    return query_records

def parse_repeat_out(repeatOut, category):
    query_records = {}
    query_lens = {}
    all_types = set()
    with open(repeatOut, 'r') as f_r:
        for line in f_r.readlines():
            if "SW" in line or "score" in line:
                continue
            info_parts = re.split(r'\s+', line.strip())
            if info_parts[0] == '':
                continue
            TE_type = info_parts[10]
            filters = ['Simple_repeat', 'Low_complexity', 'Satellite', 'Satellite/centr']
            if TE_type in filters:
                continue
            all_types.add(TE_type)

            if category is not None:
                if category == 'non-LTR':
                    if 'LINE' not in TE_type and 'SINE' not in TE_type:
                        continue
                elif category == 'Helitron':
                    if 'Helitron' not in TE_type:
                        continue
                elif category == 'LTR':
                    if 'LTR' not in TE_type:
                        continue
                elif category == 'DNA':
                    if 'DNA' not in TE_type:
                        continue

            query_name = info_parts[9].split('#')[0]
            subject_name = info_parts[4]
            strand = info_parts[8]
            if strand == "+":
                q_start = int(info_parts[11])
                q_end = int(info_parts[12])
                q_remain = int(info_parts[13].replace('(', '').replace(')', ''))
            else:
                q_start = int(info_parts[13])
                q_end = int(info_parts[12])
                q_remain = int(info_parts[11].replace('(', '').replace(')', ''))
            query_len = q_end + q_remain
            s_start = int(info_parts[5])
            s_end = int(info_parts[6])
            query_lens[query_name] = query_len
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]
            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    print(all_types)
    return query_records, query_lens

def parse_repeat_gff(repeatGff):
    query_records = {}
    with open(repeatGff, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            info_parts = line.split('\t')
            print(info_parts)
            query_name = info_parts[0].split('#')[0]
            info_list = info_parts[8].strip("\n").split(" ")
            subject_name = info_list[1].strip("Motif:").strip('"')
            q_start = int(info_list[2])
            q_end = int(info_list[3])
            s_start = int(info_parts[3])
            s_end = int(info_parts[4])
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    return query_records

def store_fasta(contigs, file_path):
    with open(file_path, 'w') as f_save:
        for name in contigs.keys():
            seq = contigs[name]
            f_save.write('>'+name+'\n'+seq+'\n')
    f_save.close()

def read_fasta(fasta_path):
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:].split(" ")[0].split('\t')[0]
                    contigseq = ''
                else:
                    contigseq += line.strip().upper()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
        rf.close()
    return contignames, contigs

def get_full_length_copies_from_blastn_v2(TE_lib, reference, alignment, tmp_output_dir, threads, type_out,
                                    full_length_threshold, category):
    ref_names, ref_contigs = read_fasta(reference)
    query_lens = {}
    if TE_lib is not None:
        query_names, query_contigs = read_fasta(TE_lib)
        for name in query_names:
            query_lens[name.split('#')[0]] = len(query_contigs[name])

    query_records = {}
    if type_out == "BLAST":
        query_records = parse_blastn_out(alignment)
    elif type_out == "RM":
        query_records, query_lens = parse_repeat_out(alignment, category)
    elif type_out == "GFF":
        query_records = parse_repeat_gff(alignment)

    all_query_copies = {}
    full_length_copies = {}
    copies_direct = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]
        if query_name not in query_lens:
            continue
        query_len = query_lens[query_name]
        skip_gap = query_len * (1 - 0.95)

        if str(query_name).__contains__('Helitron'):
            flanking_len = 5
        else:
            flanking_len = 50

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

            # cluster all closed fragments, split forward and reverse records
            forward_pos = []
            reverse_pos = []
            for pos_item in subject_pos:
                if pos_item[2] > pos_item[3]:
                    reverse_pos.append(pos_item)
                else:
                    forward_pos.append(pos_item)
            forward_pos.sort(key=lambda x: (x[2], x[3]))
            reverse_pos.sort(key=lambda x: (-x[2], -x[3]))

            clusters = {}
            cluster_index = 0
            for k, frag in enumerate(forward_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (cur_subject_start - prev_subject_end < skip_gap and cur_query_end > prev_query_end):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            cluster_index += 1
            for k, frag in enumerate(reverse_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (prev_subject_end - cur_subject_start < skip_gap and cur_query_end > prev_query_end):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            for cluster_index in clusters.keys():
                cur_cluster = clusters[cluster_index]
                cur_cluster.sort(key=lambda x: (x[0], x[1]))

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    prev_frag = cur_cluster[i]
                    if visited_frag.__contains__(prev_frag):
                        continue
                    prev_query_start = prev_frag[0]
                    prev_query_end = prev_frag[1]
                    prev_subject_start = prev_frag[2]
                    prev_subject_end = prev_frag[3]
                    prev_query_seq = (min(prev_query_start, prev_query_end), max(prev_query_start, prev_query_end))
                    prev_subject_seq = (
                        min(prev_subject_start, prev_subject_end), max(prev_subject_start, prev_subject_end))
                    prev_query_len = abs(prev_query_end - prev_query_start)
                    prev_subject_len = abs(prev_subject_end - prev_subject_start)
                    cur_longest_query_len = prev_query_len

                    cur_extend_num = 0
                    visited_frag[prev_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        cur_frag = cur_cluster[j]
                        if visited_frag.__contains__(cur_frag):
                            continue
                        cur_query_start = cur_frag[0]
                        cur_query_end = cur_frag[1]
                        cur_subject_start = cur_frag[2]
                        cur_subject_end = cur_frag[3]
                        cur_query_seq = (min(cur_query_start, cur_query_end), max(cur_query_start, cur_query_end))
                        cur_subject_seq = (min(cur_subject_start, cur_subject_end), max(cur_subject_start, cur_subject_end))

                        # could extend
                        # extend right
                        if cur_query_end > prev_query_end:
                            # judge subject direction
                            if prev_subject_start < prev_subject_end and cur_subject_start < cur_subject_end:
                                # +
                                if cur_subject_end > prev_subject_end:
                                    # forward extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and cur_subject_start - prev_subject_end < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                            elif prev_subject_start > prev_subject_end and cur_subject_start > cur_subject_end:
                                # reverse
                                if cur_subject_end < prev_subject_end:
                                    # reverse extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and prev_subject_end - cur_subject_start < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append(
                            (prev_query_start, prev_query_end, cur_longest_query_len, prev_subject_start,
                             prev_subject_end, abs(prev_subject_end - prev_subject_start), subject_name,
                             cur_extend_num))

        # To determine whether each copy has a coverage exceeding the full_length_threshold with respect
        # to the consensus sequence, retaining full-length copies.
        full_length_query_copies = {}
        full_length_flank_query_copies = {}
        query_copies = {}
        orig_query_len = query_lens[query_name]
        for repeat in longest_queries:
            # Subject
            subject_name = repeat[6]
            subject_chr_start = 0

            if repeat[3] > repeat[4]:
                direct = '-'
                old_subject_start_pos = repeat[4] - 1
                old_subject_end_pos = repeat[3]
            else:
                direct = '+'
                old_subject_start_pos = repeat[3] - 1
                old_subject_end_pos = repeat[4]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos

            subject_pos = subject_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
            subject_seq = ref_contigs[subject_name][subject_start_pos: subject_end_pos]

            flank_subject_seq = ref_contigs[subject_name][
                                subject_start_pos - flanking_len: subject_end_pos + flanking_len]
            copies_direct[subject_pos] = direct
            cur_query_len = repeat[2]
            coverage = float(cur_query_len) / orig_query_len
            if coverage >= full_length_threshold:
                full_length_query_copies[subject_pos] = subject_seq
                full_length_flank_query_copies[subject_pos] = flank_subject_seq
            query_copies[subject_pos] = (subject_name, subject_start_pos, subject_end_pos, coverage)
        full_length_copies[query_name] = full_length_query_copies
        all_query_copies[query_name] = query_copies

    # The candidate full-length copies and the consensus are then clustered using cd-hit-est,
    # retaining copies that belong to the same cluster as the consensus.
    split_files = []
    cluster_dir = tmp_output_dir + '/cluster'
    if os.path.exists(cluster_dir):
        shutil.rmtree(cluster_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    for query_name in full_length_copies.keys():
        query_copies = full_length_copies[query_name]
        fc_path = cluster_dir + '/' + query_name + '.fa'
        store_fasta(query_copies, fc_path)
        split_files.append((query_name, query_copies))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ref_index, cur_file in enumerate(split_files):
        query_name = cur_file[0]
        query_copies = cur_file[1]
        job = ex.submit(get_structure_info_v2, query_name, query_copies)
        jobs.append(job)
    ex.shutdown(wait=True)

    full_length_annotations = {}
    for job in as_completed(jobs):
        annotations = job.result()
        full_length_annotations.update(annotations)

    if os.path.exists(cluster_dir):
        shutil.rmtree(cluster_dir)
    return full_length_annotations, copies_direct, all_query_copies

def get_structure_info_v2(query_name, query_copies):
    annotations = {}
    if not annotations.__contains__(query_name):
        annotations[query_name] = []
    annotation_list = annotations[query_name]
    for copy_name in query_copies.keys():
        annotation_list.append((copy_name, ''))
    return annotations

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='Get full length TEs from alignments')
    parser.add_argument('-a', '--alignment', metavar='alignment', help='The path of alignments.')
    parser.add_argument('-l', '--lib', metavar='TE library', help='The path of TE library.')
    parser.add_argument('-r', '--ref', metavar='reference', help='The path of reference.')
    parser.add_argument('-o', '--tmp_output_dir', metavar='tmp_output_dir', help='The path of tmp_output_dir.')
    parser.add_argument('-f', '--full_length_threshold', metavar='full_length_threshold', help='The full length threshold.')
    parser.add_argument('-c', '--category', metavar='category', help='The category of TEs.')
    parser.add_argument('-s', '--species', metavar='species', help='The species.')

    args = parser.parse_args()

    alignment = args.alignment
    TE_lib = args.lib
    reference = args.ref
    tmp_output_dir = args.tmp_output_dir
    full_length_threshold = float(args.full_length_threshold)
    category = args.category
    species = args.species

    lines = generate_full_length_out_v2(alignment, TE_lib, reference, tmp_output_dir, full_length_threshold, category, type_out="RM")
    # print(lines)

    fl_bed = os.path.join(tmp_output_dir, species + '_full_length.bed')
    with open(fl_bed, 'w') as f_save:
        for line in lines:
            f_save.write(line[1] + '\t' + str(line[2]) + '\t' + str(line[3])+ '\t' + line[0] + '\n')
