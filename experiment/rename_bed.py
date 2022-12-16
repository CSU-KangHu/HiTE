import os

if __name__ == '__main__':
    name_dict = {'Chr1': 'NC_003070.9', 'Chr2': 'NC_003071.7', 'Chr3': 'NC_003074.8', 'Chr4': 'NC_003075.7', 'Chr5': 'NC_003076.8', 'ChrM': 'NC_037304.1', 'ChrC': 'NC_000932.1'}
    bed_dir = 'H:/WebTE基因文件/拟南芥/Arabidopsis_thaliana (thale cress) 3702/athaliana.ONT.50x_12345.rmet.bedgraph'
    new_bed_dir = 'H:/WebTE基因文件/拟南芥/Arabidopsis_thaliana (thale cress) 3702/athaliana.ONT.50x_12345.rmet.bedgraph.processed'
    if not os.path.exists(new_bed_dir):
        os.makedirs(new_bed_dir)

    for name in os.listdir(bed_dir):
        bed_path = bed_dir + '/' + name
        new_bed_path = new_bed_dir + '/' + name
        lines = []
        with open(bed_path, 'r') as f_r:
            for line in f_r:
                parts = line.split('\t')
                chr_name = parts[0]
                new_line = name_dict[chr_name] + '\t' + parts[1] + '\t' + parts[2] + '\t' + parts[3]
                lines.append(new_line)
        with open(new_bed_path, 'w') as f_save:
            for line in lines:
                f_save.write(line)


    # bed_path = ref_dir + '/GCF_000002035.6_GRCz11_genomic.fna'
    #
    # with open(rename_ref_path, 'w') as f_save:
    #     for name in new_contigs.keys():
    #         f_save.write('>'+name+'\n'+new_contigs[name]+'\n')