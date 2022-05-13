import argparse
import os
import random
import time

from Util import Logger, read_fasta, store_fasta

# number of sequences of different classes in model library
model_classes = {'Simple_repeat': 5, 'LTR': 10, 'DNA': 10, 'LINE': 5, 'RC/Helitron': 1}
genome_size = 10 * 1000 * 1000 # 10Mbp
min_repeat_num = 2
max_repeat_num = 30
curated_library = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/dmel_curated.fasta'
output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output'

def generate_random_genome(genome_size, genome_file):
    genome = ''
    for j in range(genome_size):
        random_base = random.choice(['A', 'T', 'C', 'G'])
        genome += random_base
    with open(genome_file, 'w') as f_save:
        f_save.write('>ref'+'\n'+genome+'\n')

def genome_simulator():
    ## step01: generate random genome sequence
    genome_file = output_dir + '/genome.fa'
    generate_random_genome(genome_size, genome_file)

    # step02: construct model library
    curated_contignames, curated_contigs = read_fasta(curated_library)
    # extract from curated library
    model_lib_path = output_dir + '/model_lib.fa'
    model_lib = {}
    for class_name in model_classes.keys():
        seq_num = model_classes[class_name]
        for name in curated_contignames:
            real_class_name = name.split('#')[1]
            if real_class_name.__contains__(class_name):
                seq = curated_contigs[name]
                seq_num -= 1
                model_lib[name] = seq
            if seq_num <= 0:
                break

    with open(model_lib_path, 'w') as f_save:
        for name in model_lib.keys():
            class_name = name.split('#')[1]
            # exclude Simple_repeat since we do not evaluate tandem repeat
            if class_name != 'Simple_repeat':
                seq = model_lib[name]
                f_save.write('>' + name + '\n' + seq + '\n')

    # step03: random insert model library into genome
    genome_contignames, genome_contigs = read_fasta(genome_file)
    genome_str = genome_contigs[genome_contignames[0]]
    print('old genome length: %d' % len(genome_str))

    for name in model_lib.keys():
        seq = model_lib[name]
        for j in range(min_repeat_num, max_repeat_num):
            insert_pos = random.choice(range(len(genome_str)))
            genome_str = genome_str[:insert_pos] + seq + genome_str[insert_pos:]
    print('new genome length: %d' % len(genome_str))

    model_genome_file = output_dir + '/genome_model.fa'
    with open(model_genome_file, 'w') as f_save:
        f_save.write('>ref' + '\n' + genome_str + '\n')

if __name__ == '__main__':
    log = Logger('GenomeSimulator.log', level='debug')

    genome_simulator()

