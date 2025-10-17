#-- coding: UTF-8 --
import argparse
import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from src.Util import read_fasta, store_fasta, Logger, split_dict_into_blocks

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run split genome chunks...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='Input genome assembly path')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')

    args = parser.parse_args()
    reference = args.g
    tmp_output_dir = args.tmp_output_dir

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    # log = Logger(tmp_output_dir + '/split.log', level='debug')

    # When passing the entire genome, the blastn alignment for obtaining TE copies consumes a large amount of memory. Therefore:
    # We divide the genome into 'threads' blocks based on chromosomes, with 'total/threads' sequences in each block.
    # Each block is of equal size and stored in a directory. Then, we pass the path as a parameter to the TIR and Helitron modules.
    # At the same time, only a portion of chromosomes is processed by one process, which does not affect the final results and also reduces memory usage.
    ref_names, ref_contigs = read_fasta(reference)
    split_ref_dir = tmp_output_dir + '/ref_chr'
    os.system('rm -rf ' + split_ref_dir)
    if not os.path.exists(split_ref_dir):
        os.makedirs(split_ref_dir)
    ref_blocks = split_dict_into_blocks(ref_contigs, 100)
    for i, block in enumerate(ref_blocks):
        new_dir = split_ref_dir + '/' + str(i)
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        chr_path = new_dir + '/ref_block_' + str(i) + '.fa'
        store_fasta(block, chr_path)
        os.system('makeblastdb -in ' + chr_path + ' -dbtype nucl' + ' > /dev/null 2>&1')





