import os
import sys


from Util import get_copies, flanking_copies, get_TSD, store_copies, multi_process_align, read_fasta, store_fasta, \
    store_copies_v1

if __name__ == '__main__':
    #我们获得的valid_tir.fa.itr文件太大，有900多M。其中有大量的冗余重复，我们根据query name取最长的一条记录即可。
    # 获取所有序列的拷贝
    threads = 48
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    reference = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    repbase_path = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase/rice/tir.repbase.ref'
    blastnResults_path = tmp_output_dir + '/test/tir.repbase.out'
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    temp_dir = tmp_output_dir + '/repbase_blast'
    multi_process_align(repbase_path, reference, blastnResults_path, blast_program_dir, temp_dir, threads)
    all_copies = get_copies(blastnResults_path, repbase_path, reference, threads=threads)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 20
    all_copies = flanking_copies(all_copies, repbase_path, reference, flanking_len, copy_num=10)
    copy_info_path = tmp_output_dir + '/test/tir.repbase.copies.info'
    store_copies_v1(all_copies, copy_info_path)

    repbase_names, repbase_contigs = read_fasta(repbase_path)

    #输出无拷贝和单拷贝序列个数
    delete_repbase_names = set()
    no_copy_num = 0
    single_copy_num = 0
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        if len(copies) == 1:
            single_copy_num += 1
            delete_repbase_names.add(query_name)
        elif len(copies) == 0:
            no_copy_num += 1
            delete_repbase_names.add(query_name)
    print('no_copy_num: '+str(no_copy_num)+', single_copy_num: '+str(single_copy_num))

    print('original repbase contig size: '+str(len(repbase_contigs)))

    for name in delete_repbase_names:
        del repbase_contigs[name]
    print('now repbase contig size: ' + str(len(repbase_contigs)))
    multi_copy_repbase_path = tmp_output_dir + '/test/tir.repbase.mc.ref'
    store_fasta(repbase_contigs, multi_copy_repbase_path)

    #找到not found的repbase序列，然后按照拷贝数排序输出
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/get_family_summary_test'
    present_file = tmp_dir + '/present.all.families'
    present_names = set()
    with open(present_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').strip()
            present_names.add(line)

    repbase_names = set(repbase_contigs.keys())
    diff_set = repbase_names.difference(present_names)
    print('not found repbase size: ' + str(len(diff_set)))

    # 取not found repbase序列的拷贝
    for name in present_names:
        if repbase_contigs.__contains__(name):
            del repbase_contigs[name]
    not_found_repbase_path = tmp_output_dir + '/test/tir.repbase.nf.ref'
    store_fasta(repbase_contigs, not_found_repbase_path)
    blastnResults_path = tmp_output_dir + '/test/tir.repbase.nf.out'
    multi_process_align(not_found_repbase_path, reference, blastnResults_path, blast_program_dir, temp_dir, threads)
    all_copies = get_copies(blastnResults_path, not_found_repbase_path, reference, threads=threads)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 20
    all_copies = flanking_copies(all_copies, not_found_repbase_path, reference, flanking_len, copy_num=10)

    # 判断copy的TSD信息
    # tsd_info = {query_name: {copy1: tsd+','+seq}, {copy2: tsd+','+seq}, {total_copy_num:}, {tsd_copy_num:}}
    tsd_info = get_TSD(all_copies, flanking_len)

    copies_list = []
    for query_name in tsd_info.keys():
        info = tsd_info[query_name]
        total_copy_num = info['total_copy_num']
        total_copy_len = info['total_copy_len']
        copies_list.append((query_name, total_copy_num, total_copy_len))
    copies_list.sort(key=lambda x: -x[2])
    print(copies_list)

    copy_info_path = tmp_output_dir + '/test/tir.repbase.nf.copies.info'
    store_copies(tsd_info, copy_info_path)


    # # 获得拷贝的数量
    # above_num = 2
    # count = 0
    # copynum = {}
    # for query_name in all_copies.keys():
    #     copies = all_copies[query_name]
    #     copy_num = len(copies)
    #     if not copynum.__contains__(copy_num):
    #         copynum[copy_num] = 0
    #     copynum[copy_num] = copynum[copy_num] + 1
    #
    #
    # for copy_num in copynum.keys():
    #     if copy_num >= above_num:
    #         count += copynum[copy_num]
    #
    # copynum = list(copynum.items())
    # copynum.sort(key=lambda x: -x[1])
    # print(copynum)
    # print(count)
    #
    # #分析单拷贝的序列，是属于什么类别，有多少个
    # single_copy = {}
    # for query_name in all_copies.keys():
    #     copies = all_copies[query_name]
    #     copy_num = len(copies)
    #     if copy_num == 1:
    #         single_copy[query_name] = 1
    #
    # repbase_path = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase/oryrep.ref'
    # class_names = {}
    # with open(repbase_path, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('>'):
    #             line = line[1:]
    #             parts = line.split('\t')
    #             query_name = parts[0]
    #             class_name = parts[1]
    #             class_names[query_name] = class_name
    #
    # count_summary = {}
    # for query_name in single_copy.keys():
    #     class_name = class_names[query_name]
    #     if not count_summary.__contains__(class_name):
    #         count_summary[class_name] = 1
    #     else:
    #         count_summary[class_name] = count_summary[class_name] + 1
    # print(count_summary)

    # #统计Repbase中超过10K的DNA转座子
    # long_TE = {}
    # contignames, contigs = read_fasta(repbase_path)
    # for name in contignames:
    #     seq = contigs[name]
    #     if len(seq) > 10000:
    #         class_name = class_names[name]
    #         if not long_TE.__contains__(class_name):
    #             long_TE[class_name] = 0
    #         long_TE[class_name] = long_TE[class_name] + 1
    # print(long_TE)