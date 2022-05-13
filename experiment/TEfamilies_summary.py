from Util import read_fasta

def print_class_num(lib_path):
    contignames, contigs = read_fasta(lib_path)
    class_names = {}
    for name in contignames:
        class_name = name.split('#')[1]
        if not class_names.__contains__(class_name):
            class_names[class_name] = 0
        num = class_names[class_name]
        class_names[class_name] = num + 1

    major_names = {}
    for name in class_names.keys():
        if name.__contains__('DNA'):
            if not major_names.__contains__('DNA'):
                major_names['DNA'] = 0
            num = major_names['DNA']
            major_names['DNA'] = num + class_names[name]
        elif name.__contains__('Helitron'):
            if not major_names.__contains__('Helitron'):
                major_names['Helitron'] = 0
            num = major_names['Helitron']
            major_names['Helitron'] = num + class_names[name]
        elif name.__contains__('LTR'):
            if not major_names.__contains__('LTR'):
                major_names['LTR'] = 0
            num = major_names['LTR']
            major_names['LTR'] = num + class_names[name]
        elif name.__contains__('LINE'):
            if not major_names.__contains__('LINE'):
                major_names['LINE'] = 0
            num = major_names['LINE']
            major_names['LINE'] = num + class_names[name]
        elif name.__contains__('SINE'):
            if not major_names.__contains__('SINE'):
                major_names['SINE'] = 0
            num = major_names['SINE']
            major_names['SINE'] = num + class_names[name]
        elif name.__contains__('Unknown'):
            if not major_names.__contains__('Unknown'):
                major_names['Unknown'] = 0
            num = major_names['Unknown']
            major_names['Unknown'] = num + class_names[name]
    print(major_names)

if __name__ == '__main__':
    krf_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib'
    #krf_libs = ['/dmel/family_dmel.fasta', '/oryza_sative/family_oryza_sativa.fasta', '/drerio/family_drerio.fasta']
    krf_libs = ['/oryza_sative/CRD.2022-05-11.11-45-41/repeats.merge.consensus.fa.final.classified', ]
    print('krf library output:')
    for lib in krf_libs:
        lib_path = krf_dir + lib
        print_class_num(lib_path)

    rm2_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/rm2_lib'
    rm2_libs = ['/Dmel-families.fa', '/rice-families.fa', '/danRer10-families.fa']
    print('rm2 library output:')
    for lib in rm2_libs:
        lib_path = rm2_dir + lib
        print_class_num(lib_path)

    curated_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib'
    curated_libs = ['/dmel_curated.fasta', '/O.sativa_curated.fasta', '/D.rerio_curated.fasta']
    print('curated library output:')
    for lib in curated_libs:
        lib_path = curated_dir + lib
        print_class_num(lib_path)