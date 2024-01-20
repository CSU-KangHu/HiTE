import argparse
import os
import random
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from configs import config
from utils.data_util import read_fasta_v1, read_fasta, store_fasta, extract_non_autonomous


def main():
    # 1.parse args
    describe_info = '########################## NeuralTE-extract_species, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data_path', metavar='data_path', help='Input Repbase data path')
    parser.add_argument('--out_dir', metavar='out_dir', help='Output directory for training set (non-xxx species) and testing set (xxx species)')
    parser.add_argument('--species', metavar='species', help='Input species to be testing set')

    args = parser.parse_args()

    data_path = args.data_path
    out_dir = args.out_dir
    species = args.species

    data_path = os.path.realpath(data_path)
    out_dir = os.path.realpath(out_dir)

    work_dir = out_dir
    train = work_dir + '/train.ref'
    test = work_dir + '/test.ref'
    names, contigs = read_fasta_v1(data_path)
    train_contigs = {}
    test_contigs = {}
    for name in names:
        species_name = name.split('\t')[2]
        label = name.split('\t')[1]
        if species_name == species and label != 'Unknown':
            test_contigs[name] = contigs[name]
        elif not species_name.__contains__(species.split(' ')[0]):
            train_contigs[name] = contigs[name]
    store_fasta(train_contigs, train)
    store_fasta(test_contigs, test)

if __name__ == '__main__':
    main()