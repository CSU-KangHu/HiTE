import os
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

import argparse

from configs import config
from utils.data_util import read_fasta, read_fasta_v1, store_fasta


def main():
    # 1.parse args
    describe_info = '########################## NeuralTE-reFormat_RM2, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('-i', metavar='input_fasta', help='Input the file path to be reformatted, fasta format')
    parser.add_argument('-o', metavar='output_fasta', help='Output the reformatted file, fasta format')

    args = parser.parse_args()

    data_path = args.i
    output_path = args.o
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    names, contigs = read_fasta_v1(data_path)
    intact_LTRs = {}
    intact_LTRs_set = set()
    for name in names:
        if name.__contains__('Type=LTR'):
            family_num = int(name.split('#')[0].split('ltr-1_')[1].split('-')[1])
            if not intact_LTRs.__contains__(family_num):
                intact_LTRs[family_num] = []
            elements = intact_LTRs[family_num]
            elements.append(name)
            intact_LTRs_set.add(name)
    for name in names:
        if name.__contains__('Type=INT'):
            family_num = int(name.split('#')[0].split('ltr-1_')[1].split('-')[1]) - 1
            if not intact_LTRs.__contains__(family_num):
                intact_LTRs[family_num] = []
            elements = intact_LTRs[family_num]
            elements.append(name)
            intact_LTRs_set.add(name)
    new_contigs = {}
    for name in names:
        if name not in intact_LTRs_set:
            new_contigs[name] = contigs[name]
    for family_num in intact_LTRs.keys():
        for i, name in enumerate(intact_LTRs[family_num]):
            if i == 0:
                type = 'LTR'
            else:
                type = 'I'
            new_name = 'RM2_intactLTR_' + str(family_num) + '-' + type
            new_contigs[new_name] = contigs[name]
    store_fasta(new_contigs, output_path)

if __name__ == '__main__':
    main()