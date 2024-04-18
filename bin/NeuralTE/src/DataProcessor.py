# This program loads data and determines whether data needs respective preprocessing
# based on which features are used in the configuration file.
import os
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from configs import config
from concurrent.futures import ProcessPoolExecutor, as_completed

from utils.data_util import read_fasta_v1, generate_TSD_info, generate_domain_info, generate_terminal_info, \
    load_repbase_with_TSD, generate_feature_mats, read_fasta, connect_LTR, store_fasta, generate_minority_info


class DataProcessor:
    def __init__(self):
        self.project_dir = config.project_dir
        self.tool_dir = self.project_dir + '/tools'
        self.work_dir = config.work_dir
        self.threads = config.threads
        self.minority_labels_class = config.minority_labels_class
        self.all_wicker_class = config.all_wicker_class
        self.use_terminal = config.use_terminal
        self.use_TSD = config.use_TSD
        self.use_domain = config.use_domain
        self.use_minority = config.use_minority
        self.use_ends = config.use_ends
        self.is_train = config.is_train
        self.is_predict = config.is_predict

        self.ex = ProcessPoolExecutor(self.threads)


    def load_data(self, internal_kmer_sizes, terminal_kmer_sizes, data_path):
        # Copy input files to the working directory
        if not os.path.exists(data_path):
            print('Input file not exist: ' + data_path)
            exit(-1)
        os.makedirs(config.work_dir, exist_ok=True)
        os.system('cp ' + data_path + ' ' + config.work_dir)
        genome_info_path = config.work_dir + '/genome.info'
        data_path = config.work_dir + '/' + os.path.basename(data_path)
        domain_train_path = data_path + '.domain'

        minority_temp = config.work_dir + '/minority'
        if not os.path.exists(minority_temp):
            os.makedirs(minority_temp)
        minority_train_path = minority_temp + '/train.minority.ref'
        minority_out = minority_temp + '/train.minority.out'

        data_path = self.preprocess_data(data_path, domain_train_path, minority_train_path, minority_out, self.work_dir,
                                                    self.project_dir, self.tool_dir, self.threads,
                                                    self.use_TSD, self.use_domain, self.use_minority, self.use_terminal,
                                                    self.is_predict, self.is_train, genome_info_path)

        X, Y, seq_names = load_repbase_with_TSD(data_path, domain_train_path, minority_train_path, minority_out,
                                                    self.all_wicker_class, self.project_dir+'/data/TEClasses.tsv')

        X, Y = generate_feature_mats(X, Y, seq_names, self.minority_labels_class, self.all_wicker_class,
                                                    internal_kmer_sizes, terminal_kmer_sizes, self.ex)

        # Reshape data into the format accepted by the model
        X = X.reshape(X.shape[0], config.X_feature_len, 1)
        X = X.astype('float64')
        return X, Y, seq_names, data_path

    def preprocess_data(self, data, domain_train_path, minority_train_path, minority_out, work_dir, project_dir,
                        tool_dir, threads, use_TSD, use_domain, use_minority, use_terminal, is_predict, is_train, genome_info_path):
        # Delete previous run's retained results
        SegLTR2intactLTRMap = config.work_dir + '/segLTR2intactLTR.map'
        os.system('rm -f ' + SegLTR2intactLTRMap)

        if is_predict:
            # Format input files to Repbase format
            names, contigs = read_fasta(data)
            os.makedirs(os.path.dirname(data), exist_ok=True)
            with open(data, 'w') as f_save:
                for name in names:
                    seq = contigs[name]
                    name = name.split('/')[0].split('#')[0]
                    new_name = name + '\tUnknown\tUnknown'
                    f_save.write('>' + new_name + '\n' + seq + '\n')

            # Concatenate LTR sequences with LTR internal sequences from the input TE library to create a complete LTR sequence
            data, repbase_labels = connect_LTR(data)
            # Remove domain files to ensure regeneration every time
            os.system('rm -f ' + domain_train_path)
            # set keep_raw = 1 to ensure not lose any sequence
            config.keep_raw = 1

        names, contigs = read_fasta_v1(data)
        # Convert Repbase labels to Wicker format
        converted_contigs = {}
        unconverted_contigs = {}
        for name in names:
            parts = name.split('\t')
            if len(parts) >= 3:
                label = parts[1]
                if config.all_wicker_class.__contains__(label):
                    new_label = label
                elif config.Repbase_wicker_labels.__contains__(label):
                    new_label = config.Repbase_wicker_labels[label]
                else:
                    new_label = None
                if new_label is not None:
                    new_name = '\t'.join([parts[0], new_label] + parts[2:])
                    converted_contigs[new_name] = contigs[name]
                else:
                    unconverted_contigs[name] = contigs[name]
            else:
                unconverted_contigs[name] = contigs[name]
        store_fasta(converted_contigs, data)
        if len(unconverted_contigs) > 0:
            unconverted_data = config.work_dir + '/unconverted_TE.fa'
            store_fasta(unconverted_contigs, unconverted_data)
            print('Warning: The input TE library contains unknown superfamily labels, total size = ' + str(len(unconverted_contigs)) + ', which saved at ' + os.path.realpath(unconverted_data))

        # Retrieve a record and check if it contains the corresponding information; if not, it needs to be processed accordingly
        if len(names) > 0:
            cur_name = names[0]
            if use_TSD != 0 and ('TSD:' not in cur_name or 'TSD_len:' not in cur_name):
                # Use TSD feature, but if TSD information is missing in the dataset, it needs to be regenerated
                # Ensure modifying the 'data/ncbi_ref.info' file with the correct Genome Path
                is_expanded = 0 # Set is_expanded=1 only when balancing the Repbase dataset; otherwise, set is_expanded=0
                keep_raw = config.keep_raw # Retain original sequences; set to 0 to retain only sequences with TSD, useful when training the TSD model.
                print('keep_raw:' + str(keep_raw))
                data = generate_TSD_info(data, genome_info_path, work_dir, is_expanded, keep_raw, threads)
            if use_domain != 0 and not os.path.exists(domain_train_path):
                # Use domain feature, but if the domain information file is missing, it needs to be regenerated through comparison.
                generate_domain_info(data, project_dir+'/data/RepeatPeps.lib', work_dir, threads)
            if use_minority != 0:
                generate_minority_info(data, minority_train_path, minority_out, threads, is_train)
            if use_terminal != 0 and ('LTR:' not in cur_name or 'TIR:' not in cur_name):
                data = generate_terminal_info(data, work_dir, tool_dir, threads)
        return data
