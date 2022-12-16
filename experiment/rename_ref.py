from Util import read_fasta

if __name__ == '__main__':
    ref_dir = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/drerio'
    ref_path = ref_dir + '/GCF_000002035.6_GRCz11_genomic.fna'
    ref_dir1 = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/drerio'
    rename_ref_path = ref_dir1 + '/GCF_000002035.6_GRCz11_genomic.fna'
    ref_contigNames, ref_contigs = read_fasta(ref_path)
    new_contigs = {}
    for id, name in enumerate(ref_contigNames):
        ref_seq = ref_contigs[name]
        new_contigs['chr_'+str(id)] = ref_seq

    with open(rename_ref_path, 'w') as f_save:
        for name in new_contigs.keys():
            f_save.write('>'+name+'\n'+new_contigs[name]+'\n')