import os

from utils.data_util import find_files_recursively, read_fasta_v1, store_fasta

if __name__ == '__main__':
    # 从Repbase目录下获取某些物种的金标准数据集
    work_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/RepBase28.06.fasta'
    species_name = 'Thalassiosira pseudonana'
    tmp_out_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/Repbase_28.06/Thalassiosira_pseudonana'

    file_extension = '.ref'
    ref_files = find_files_recursively(work_dir, file_extension)
    print(ref_files)

    tags = set()
    repbase_names = []
    repbase_contigs = {}
    for cur_file in ref_files:
        ref_names, ref_contigs = read_fasta_v1(cur_file)
        for name in ref_names:
            if species_name in name:
                repbase_names.append(name)
                repbase_contigs[name] = ref_contigs[name]
                tag = name.split('\t')[1]
                tags.add(tag)
    print(tags)
    print(len(tags))

    ltr_tags = ['Gypsy', 'Copia', 'LTR Retrotransposon', 'BEL', 'LTR', 'Endogenous Retrovirus', 'Caulimoviridae', 'Long terminal repeat']
    tir_tags = ['Mariner/Tc1', 'DNA transposon', 'EnSpm/CACTA', 'MuDR', 'hAT', 'Harbinger', 'Transib', 'piggyBac', 'P',
                'DNA', 'Sola2', 'Kolobok', ]
    helitron_tags = ['Helitron', 'MINIME_DN']
    non_ltr_tags = ['L1', 'SINE2/tRNA', 'Non-LTR Retrotransposon', 'SINE', 'R1', 'Jockey', 'CR1', 'R2', 'RTEX', 'Hero',
                    'RTE', 'tRNA']

    if not os.path.exists(tmp_out_dir):
        os.makedirs(tmp_out_dir)
    ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.repbase.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.repbase.ref'
    all_repbase_path = tmp_out_dir + '/all.repbase.ref'

    all_contigs = {}
    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    other_tags = set()
    for name in repbase_names:
        tag = name.split('\t')[1]
        if tag in ltr_tags:
            ltr_contigs[name] = repbase_contigs[name]
        elif tag in tir_tags:
            tir_contigs[name] = repbase_contigs[name]
        elif tag in helitron_tags:
            helitron_contigs[name] = repbase_contigs[name]
        elif tag in non_ltr_tags:
            non_ltr_contigs[name] = repbase_contigs[name]
        else:
            other_tags.add(tag)
        all_contigs[name] = repbase_contigs[name]
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)
    store_fasta(all_contigs, all_repbase_path)
    print(other_tags)