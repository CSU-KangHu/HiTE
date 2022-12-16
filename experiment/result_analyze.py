from Util import read_fasta

if __name__ == '__main__':
    dir_path = '/home/hukang/KRF_output/drerio/CRD.2022-08-21.19-16-38'
    file_path = dir_path + '/longest_repeats.confident.cons.fa'
    test_path = dir_path + '/module.fa'
    names, contigs = read_fasta(file_path)
    counts = {}
    for name in names:
        seq = contigs[name]
        parts = name.split('-')
        from_type = parts[2].split('#')[0]
        if not counts.__contains__(from_type):
            counts[from_type] = 0
        ori_count = counts[from_type]
        counts[from_type] = ori_count + 1
    print(counts)

    rm2_path = '/home/hukang/KmerRepFinder_test/library/rm2_lib/danRer10-families.fa'
    names, contigs = read_fasta(rm2_path)
    counts = {}
    for name in names:
        parts = name.split('#')
        class_name = parts[1].split('/')[0]
        if not counts.__contains__(class_name):
            counts[class_name] = 0
        ori_count = counts[class_name]
        counts[class_name] = ori_count + 1
    print(counts)