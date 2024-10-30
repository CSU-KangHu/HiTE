#-- coding: UTF-8 --
import argparse
import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, convertToUpperCase_v1, multi_line, split_dict_into_blocks

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE split genome chunks...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='Input genome assembly path')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--chrom_seg_length', metavar='chrom_seg_length',
                        help='The length of genome segments')
    parser.add_argument('--chunk_size', metavar='chunk_size',
                        help='The chunk size of large genome')

    args = parser.parse_args()
    reference = args.g
    tmp_output_dir = args.tmp_output_dir
    chrom_seg_length = int(args.chrom_seg_length)
    chunk_size = float(args.chunk_size)

    chunk_size = int(chunk_size * 1024 * 1024)

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir + '/HiTE_split.log', level='debug')

    log.logger.info('Start Splitting Reference into chunks')
    reference_pre = convertToUpperCase_v1(reference)
    reference_tmp = multi_line(reference_pre, chrom_seg_length)
    cut_references = []
    cur_ref_contigs = {}
    cur_base_num = 0
    ref_index = 0
    with open(reference_tmp, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ref_name = parts[0].replace('>', '')
            start = parts[1]
            seq = parts[2]
            new_ref_name = ref_name + '$' + start
            cur_ref_contigs[new_ref_name] = seq
            cur_base_num += len(line)
            if cur_base_num >= chunk_size:
                # store references
                cur_ref_path = tmp_output_dir + '/genome.cut' + str(ref_index) + '.fa'
                store_fasta(cur_ref_contigs, cur_ref_path)
                cut_references.append(cur_ref_path)
                cur_ref_contigs = {}
                cur_base_num = 0
                ref_index += 1
        if len(cur_ref_contigs) > 0:
            cur_ref_path = tmp_output_dir + '/genome.cut' + str(ref_index) + '.fa'
            store_fasta(cur_ref_contigs, cur_ref_path)
            cut_references.append(cur_ref_path)
    f_r.close()

    # When passing the entire genome, the blastn alignment for obtaining TE copies consumes a large amount of memory. Therefore:
    # We divide the genome into 'threads' blocks based on chromosomes, with 'total/threads' sequences in each block.
    # Each block is of equal size and stored in a directory. Then, we pass the path as a parameter to the TIR and Helitron modules.
    # At the same time, only a portion of chromosomes is processed by one process, which does not affect the final results and also reduces memory usage.
    ref_names, ref_contigs = read_fasta(reference)
    split_ref_dir = tmp_output_dir + '/ref_chr'
    os.system('rm -rf ' + split_ref_dir)
    if not os.path.exists(split_ref_dir):
        os.makedirs(split_ref_dir)
    ref_blocks = split_dict_into_blocks(ref_contigs, 100, chunk_size)
    for i, block in enumerate(ref_blocks):
        chr_path = split_ref_dir + '/ref_block_' + str(i) + '.fa'
        store_fasta(block, chr_path)
        os.system('makeblastdb -in ' + chr_path + ' -dbtype nucl'+ ' > /dev/null 2>&1')





