import os
import sys



cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta_v1, store_fasta
from main import calculate_max_min

def get_score(confident_TIR):
    # (tir_len, tsd, te_len, contigs[name])
    TE_len_list = []
    tir_len_list = []
    tsd_len_list = []
    for info in confident_TIR:
        TE_len_list.append(info[2])
        tir_len_list.append(info[0])
        tsd_len_list.append(len(info[1]))

    max_TE_len = max(TE_len_list)
    min_TE_len = min(TE_len_list)

    max_tir_len = max(tir_len_list)
    min_tir_len = min(tir_len_list)

    max_tsd_len = max(tsd_len_list)
    min_tsd_len = min(tsd_len_list)

    max_info = None
    max_score = -1
    for info in confident_TIR:
        if max_TE_len == min_TE_len:
            TE_len_normal = 0
        else:
            TE_len_normal = calculate_max_min(info[2], max_TE_len, min_TE_len)
        if max_tir_len == min_tir_len:
            tir_len_normal = 0
        else:
            tir_len_normal = calculate_max_min(info[0], max_tir_len, min_tir_len)
        if max_tsd_len == min_tsd_len:
            tsd_len_normal = 0
        else:
            tsd_len_normal = calculate_max_min(len(info[1]), max_tsd_len, min_tsd_len)
        score = 0.3*TE_len_normal+0.3*tir_len_normal+0.4*tsd_len_normal
        if score > max_score:
            max_score = score
            max_info = info
    return max_info

if __name__ == '__main__':
    #我们获得的valid_tir.fa.itr文件太大，有900多M。其中有大量的冗余重复，我们根据query name取最长的一条记录即可。
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    input_path = tmp_dir + '/tir_tsd.fa.itr'
    output_path = tmp_dir + '/tir_tsd.filter_dup.fa'

    #综合考虑TSD len, tir len, TE len等多个因素
    contignames, contigs = read_fasta_v1(input_path)
    filtered_contigs = {}
    for name in contignames:
        tir_len = int(name.split('Length itr=')[1])
        parts = name.split('-C_')
        orig_query_name = parts[0]
        tsd = parts[1].split(' ')[0].split('-tsd_')[1]
        te_len = len(contigs[name])
        #te_len = int(parts[1].split('-')[1].split('_')[1])
        if not filtered_contigs.__contains__(orig_query_name):
            filtered_contigs[orig_query_name] = set()
        confident_TIR = filtered_contigs[orig_query_name]
        confident_TIR.add((tir_len, tsd, te_len, contigs[name]))

    node_index = 0
    with open(output_path, 'w') as f_save:
        for name in filtered_contigs.keys():
            confident_TIR = filtered_contigs[name]
            highest_confident_TIR = get_score(confident_TIR)
            query_name = 'N_'+str(node_index)+'-tirlen_'+str(highest_confident_TIR[0])+'-TElen_'+str(highest_confident_TIR[2])+'-tsd_'+str(highest_confident_TIR[1])
            f_save.write('>'+query_name+'\n'+highest_confident_TIR[3]+'\n')
            node_index += 1
    f_save.close()




    #挑选TElen最长的序列来代表最终的TIR序列
    # contignames, contigs = read_fasta_v1(input_path)
    # filtered_contigs = {}
    # for name in contignames:
    #     tir_len = int(name.split('Length itr=')[1])
    #     parts = name.split('-C_')
    #     orig_query_name = parts[0]
    #     tsd = parts[1].split(' ')[0].split('-tsd_')[1]
    #     #te_len = int(parts[1].split('-')[1].split('_')[1])
    #     if not filtered_contigs.__contains__(orig_query_name):
    #         te_len = len(contigs[name])
    #         filtered_contigs[orig_query_name] = (tir_len, tsd, te_len, contigs[name])
    #     else:
    #         old_seq_item = filtered_contigs[orig_query_name]
    #         new_seq = contigs[name]
    #         if len(new_seq) > old_seq_item[2]:
    #             filtered_contigs[orig_query_name] = (tir_len, tsd, len(new_seq), new_seq)
    #
    # node_index = 0
    # with open(output_path, 'w') as f_save:
    #     for name in filtered_contigs.keys():
    #         item = filtered_contigs[name]
    #         query_name = 'N_'+str(node_index)+'-tirlen_'+str(item[0])+'-TElen_'+str(item[2])+'-tsd_'+str(item[1])
    #         f_save.write('>'+query_name+'\n'+item[3]+'\n')
    #         node_index += 1
    # f_save.close()

    #335 perfect.families