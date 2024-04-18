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
    describe_info = '########################## NeuralTE-extract_non_auto_DNA, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data_path', metavar='data_path', help='Input Repbase data path')
    parser.add_argument('--out_path', metavar='out_path', help='Output non-autonomous DNA transposons')

    args = parser.parse_args()

    data_path = args.data_path
    out_path = args.out_path

    data_path = os.path.realpath(data_path)
    out_path = os.path.realpath(out_path)

    # 抽出非自治转座子
    extract_non_autonomous(data_path, out_path)

if __name__ == '__main__':
    main()