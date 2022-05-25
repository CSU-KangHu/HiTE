import argparse
import os
import random
import time

from Util import Logger, read_fasta, store_fasta

# number of sequences of different classes in model library
#model_classes = {'Simple_repeat': 5, 'LTR': 10, 'DNA': 10, 'LINE': 5, 'RC/Helitron': 1}
genome_size = 100 * 1000 * 1000 # 100Mbp
min_repeat_num = 2
max_repeat_num = 75
max_mutation_ratio = 0.1
curated_library = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/dmel_curated.fasta'
output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output_2-75'

def generate_random_genome(genome_size, genome_file):
    genome = ''
    for j in range(genome_size):
        random_base = random.choice(['A', 'T', 'C', 'G'])
        genome += random_base
    with open(genome_file, 'w') as f_save:
        f_save.write('>ref'+'\n'+genome+'\n')


def genome_simulator():
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ## step01: generate random genome sequence
    genome_file = output_dir + '/genome.fa'
    generate_random_genome(genome_size, genome_file)

    # step02: construct model library
    curated_contignames, curated_contigs = read_fasta(curated_library)

    # # extract from curated library
    model_lib_path = output_dir + '/model_lib.fa'
    with open(model_lib_path, 'w') as f_save:
        for name in curated_contignames:
            if name == 'PROTOP_B#DNA/P':
                continue
            seq = curated_contigs[name]
            f_save.write('>' + name + '\n' + seq + '\n')
    #os.system('cp ' + curated_library + ' ' + model_lib_path)
    model_names, model_lib = read_fasta(curated_library)
    # model_lib = {}
    # for class_name in model_classes.keys():
    #     seq_num = model_classes[class_name]
    #     for name in curated_contignames:
    #         real_class_name = name.split('#')[1]
    #         if real_class_name.__contains__(class_name):
    #             seq = curated_contigs[name]
    #             seq_num -= 1
    #             model_lib[name] = seq
    #         if seq_num <= 0:
    #             break
    #
    # with open(model_lib_path, 'w') as f_save:
    #     for name in model_lib.keys():
    #         class_name = name.split('#')[1]
    #         # exclude Simple_repeat since we do not evaluate tandem repeat
    #         seq = model_lib[name]
    #         f_save.write('>' + name + '\n' + seq + '\n')

    # step03: random insert model library into genome
    genome_contignames, genome_contigs = read_fasta(genome_file)
    genome_str = genome_contigs[genome_contignames[0]]
    print('old genome length: %d' % len(genome_str))
    print('total model library size = %d' %len(model_lib.keys()))
    for i, name in enumerate(model_lib.keys()):
        seq = model_lib[name]
        seq = simulate_mutation(seq)
        repeat_num = random.randint(min_repeat_num, max_repeat_num)
        print('insert %d th sequence of curated library, repeat num = %d' %(i+1, repeat_num))
        for j in range(repeat_num):
            insert_pos = random.choice(range(len(genome_str)))
            genome_str = genome_str[:insert_pos] + seq + genome_str[insert_pos:]
    print('new genome length: %d' % len(genome_str))

    model_genome_file = output_dir + '/genome_model.fa'
    with open(model_genome_file, 'w') as f_save:
        f_save.write('>ref' + '\n' + genome_str + '\n')

def random_insertion(total_insertion_len, seq):
    random_list = []
    random_len = 0
    while random_len < total_insertion_len:
        random_num = random.randint(1, total_insertion_len)
        random_len += random_num
        if random_len <= total_insertion_len:
            random_list.append(random_num)
        else:
            random_list.append(total_insertion_len - random_len + random_num)

    for insertion_size in random_list:
        insertion_str = ''
        for j in range(insertion_size):
            random_base = random.choice(['A', 'T', 'C', 'G'])
            insertion_str += random_base
        insert_pos = random.choice(range(len(seq)))
        seq = seq[:insert_pos] + insertion_str + seq[insert_pos:]
    return seq


def random_snp(total_snp_len, seq):
    seq_list = list(seq)
    seq_positions = range(len(seq_list))
    snped_pos = []
    for i in range(total_snp_len):
        snp_pos = random.choice(list(set(seq_positions)-set(snped_pos)))
        random_base = random.choice(list(set(['A', 'T', 'C', 'G']) - set(seq_list[snp_pos])))
        seq_list[snp_pos] = random_base
        snped_pos.append(snp_pos)
    return ''.join(seq_list)


def random_deletion(total_deletion_len, seq):
    del_len = 0
    seq_list = list(seq)
    del_positions = range(len(seq_list))
    deled_pos = []
    while del_len < total_deletion_len:
        del_pos = random.choice(list(set(del_positions)-set(deled_pos)))
        seq_list[del_pos] = ''
        del_len += 1
        deled_pos.append(del_pos)
    return ''.join(seq_list)

def simulate_mutation(seq):
    mutatio_ratio = random.uniform(0, max_mutation_ratio)
    mutatio_len = int(len(seq) * mutatio_ratio)
    total_insertion_len = int(mutatio_len/3)
    total_snp_len = int(mutatio_len/3)
    total_deletion_len = mutatio_len - (total_insertion_len + total_snp_len)
    #print('simulate insertion size: %d' %total_insertion_len)
    seq = random_insertion(total_insertion_len, seq)
    #print('simulate snp size: %d' % total_snp_len)
    seq = random_snp(total_snp_len, seq)
    #print('simulate deletion size: %d' % total_deletion_len)
    seq = random_deletion(total_deletion_len, seq)
    return seq

if __name__ == '__main__':
    log = Logger('GenomeSimulator.log', level='debug')

    genome_simulator()


