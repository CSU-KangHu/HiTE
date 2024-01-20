import argparse
import os
import random
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from configs import config
from utils.data_util import read_fasta_v1, read_fasta, store_fasta

def split_dataset(sequences, train_output_file, test_output_file, split_ratio=0.8):
    # Save the dataset by category, dividing data in each category into 80% for training and 20% for testing.
    sequences_types = {}
    for name in sequences.keys():
        label = name.split('\t')[1]
        if not sequences_types.__contains__(label):
            sequences_types[label] = {}
        cur_label_sequences = sequences_types[label]
        cur_label_sequences[name] = sequences[name]

    train_sequences = {}
    test_sequences = {}
    for label in sequences_types.keys():
        cur_label_sequences = sequences_types[label]

        # Randomly partition the dataset into training and testing sets.
        all_ids = list(cur_label_sequences.keys())
        random.shuffle(all_ids)
        train_ids = all_ids[:int(split_ratio * len(all_ids))]
        test_ids = all_ids[int(split_ratio * len(all_ids)):]
        for seq_id in train_ids:
            train_sequences[seq_id] = cur_label_sequences[seq_id]
        for seq_id in test_ids:
            test_sequences[seq_id] = cur_label_sequences[seq_id]

    store_fasta(train_sequences, train_output_file)
    store_fasta(test_sequences, test_output_file)


def print_dataset_info(repbase_path, type):
    repbase_names, repbase_contigs = read_fasta_v1(repbase_path)
    # Count the number of sequences and species within the dataset.
    unique_species = set()

    for name in repbase_names:
        parts = name.split('\t')
        unique_species.add(parts[2])
    if type == 'train':
        print('Train dataset sequence size: ' + str(len(repbase_names)) + ', total species num: ' + str(len(unique_species)))
    else:
        print('Test dataset sequence size: ' + str(len(repbase_names)) + ', total species num: ' + str(len(unique_species)))


def main():
    # 1.parse args
    describe_info = '########################## NeuralTE-split_train_test, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data_path', metavar='data_path', help='Input data path to be splitted')
    parser.add_argument('--out_dir', metavar='out_dir', help='Output directory')
    parser.add_argument('--ratio', metavar='ratio', help='Ratio of training set to test set')

    args = parser.parse_args()

    data_path = args.data_path
    out_dir = args.out_dir
    ratio = args.ratio

    if ratio is None:
        ratio = 0.8
    else:
        ratio = float(ratio)

    data_path = os.path.realpath(data_path)
    out_dir = os.path.realpath(out_dir)
    repbase_train_path = out_dir + '/train.ref'
    repbase_test_path = out_dir + '/test.ref'

    # Divide the dataset into training and testing sets in an 80:20 ratio.
    repbase_names, repbase_contigs = read_fasta_v1(data_path)
    split_dataset(repbase_contigs, repbase_train_path, repbase_test_path, split_ratio=ratio)
    print_dataset_info(repbase_train_path, 'train')
    print_dataset_info(repbase_test_path, 'test')


if __name__ == '__main__':
    main()