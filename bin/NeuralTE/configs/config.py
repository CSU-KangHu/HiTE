#-- coding: UTF-8 --
# config.py：
# This file stores variables and parameters defined in NeuralTE.
# Ensure you understand its purpose before making modifications;
# otherwise, it is advisable to keep the default values.

import os
from multiprocessing import cpu_count
current_folder = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.join(current_folder, "..")

# 1. Data preprocessing parameters
## Whether to use corresponding features for classification, all of which have been proven helpful for classification.
use_kmers = 1   # Whether to use k-mer feature
use_terminal = 1    # Whether to use LTR and TIR features
use_TSD = 0     # Whether to use TSD feature
use_domain = 1  # Whether to use TE domain feature
use_ends = 1    # Whether to use 5-bp ends feature


use_minority = 0 # Whether to use minority samples to correct results
is_train = 0  # Whether it is in the model training stage
keep_raw = 0  # Whether to retain the raw input sequence, 1 yes, 0 no, only save species having TSDs
only_preprocess = 0 # Whether to only perform data preprocessing
is_predict = 1  # Enable prediction mode. Setting to 0 requires the input FASTA file to be in Repbase format (seq_name\tLabel\tspecies).
is_wicker = 1   # Use Wicker classification labels. Setting to 0 will output RepeatMasker classification labels.
is_plant = 0 # Is the input genome of a plant? 0 represents non-plant, while 1 represents plant.
is_debug = 0 # Is debug mode


# 2. Program and model parameters
threads = int(cpu_count())  # Use multi-threading to load data
internal_kmer_sizes = [5]   # Size of k-mer used for converting internal sequences to k-mer frequency features
terminal_kmer_sizes = [3, 4] # Size of k-mer used for converting terminal sequences to k-mer frequency features
## CNN model parameters
cnn_num_convs = 3 # Number of CNN convolutional layers
cnn_filters_array = [16, 16, 16] # Number of filters per convolutional layer in CNN
cnn_kernel_sizes_array = [7, 7, 7] # Kernel size for each convolutional layer in CNN; for 2D convolutional layers, set as [(3, 3), ...]
cnn_dropout = 0.5 # CNN dropout threshold
## Training parameters
batch_size = 32 # Batch size for training
epochs = 50 #  Number of epochs for training
use_checkpoint = 0  # Whether to use checkpoint training; set to 1 to resume training from the parameters of the previous failed training session, avoiding training from scratch


################################################### The following parameters do not need modification ######################################################################
version_num = '1.0.0'
work_dir = project_dir + '/work' # temp work directory

non_temp_files = ['classified\.info', 'classified_TE\.fa', '.*\.domain']

# minority sample labels
#minority_labels_class = {'Crypton': 0, '5S': 1, '7SL': 2, 'Merlin': 3, 'P': 4, 'R2': 5, 'Unknown': 6}
minority_labels_class = {'Crypton': 0, '5S': 1, 'Merlin': 2, 'P': 3, 'R2': 4, 'Unknown': 5}

## Superfamily labels based on Wicker classification
all_wicker_class = {'Tc1-Mariner': 0, 'hAT': 1, 'Mutator': 2, 'Merlin': 3, 'Transib': 4, 'P': 5, 'PiggyBac': 6,
                    'PIF-Harbinger': 7, 'CACTA': 8, 'Crypton': 9, 'Helitron': 10, 'Maverick': 11, 'Copia': 12,
                    'Gypsy': 13, 'Bel-Pao': 14, 'Retrovirus': 15, 'DIRS': 16, 'Ngaro': 17, 'VIPER': 18,
                    'Penelope': 19, 'R2': 20, 'RTE': 21, 'Jockey': 22, 'L1': 23, 'I': 24, 'tRNA': 25, '7SL': 26, '5S': 27, 'Unknown': 28}
class_num = len(all_wicker_class)
inverted_all_wicker_class = {value: key for key, value in all_wicker_class.items()}
# Maximum length of TSD (Target Site Duplication)
max_tsd_length = 15
# Obtain CNN input dimensions
X_feature_len = 0
# Dimensions of TE terminal and internal sequences
if use_kmers != 0:
    for kmer_size in internal_kmer_sizes:
        X_feature_len += pow(4, kmer_size)
    if use_terminal != 0:
        for i in range(2):
            for kmer_size in terminal_kmer_sizes:
                X_feature_len += pow(4, kmer_size)
if use_TSD != 0:
    X_feature_len += max_tsd_length * 4 + 1
# if use_minority != 0:
#     X_feature_len += len(minority_labels_class)
if use_domain != 0:
    X_feature_len += len(all_wicker_class)
if use_ends != 0:
    X_feature_len += 10 * 4

## Mapping Repbase labels to Wicker labels
Repbase_wicker_labels = {'Mariner/Tc1': 'Tc1-Mariner', 'mariner/Tc1 superfamily': 'Tc1-Mariner', 'hAT': 'hAT',
                         'HAT superfamily': 'hAT', 'MuDR': 'Mutator', 'Merlin': 'Merlin', 'Transib': 'Transib',
                         'P': 'P', 'P-element': 'P', 'PiggyBac': 'PiggyBac', 'Harbinger': 'PIF-Harbinger',
                         'EnSpm/CACTA': 'CACTA', 'Crypton': 'Crypton', 'CryptonF': 'Crypton', 'CryptonS': 'Crypton',
                         'CryptonI': 'Crypton', 'CryptonV': 'Crypton', 'CryptonA': 'Crypton', 'Helitron': 'Helitron',
                         'HELITRON superfamily': 'Helitron', 'Copia': 'Copia', 'Gypsy': 'Gypsy',
                         'GYPSY superfamily': 'Gypsy', 'Gypsy retrotransposon': 'Gypsy', 'BEL': 'Bel-Pao',
                         'Caulimoviridae': 'Retrovirus',
                         'ERV1': 'Retrovirus', 'ERV2': 'Retrovirus', 'ERV3': 'Retrovirus', 'Lentivirus': 'Retrovirus',
                         'ERV4': 'Retrovirus', 'Lokiretrovirus': 'Retrovirus', 'DIRS': 'DIRS', 'Penelope': 'Penelope',
                         'Penelope/Poseidon': 'Penelope', 'Neptune': 'Penelope', 'Nematis': 'Penelope',
                         'Athena': 'Penelope', 'Coprina': 'Penelope', 'Hydra': 'Penelope', 'Naiad/Chlamys': 'Penelope',
                         'R2': 'R2', 'RTE': 'RTE', 'Jockey': 'Jockey', 'L1': 'L1', 'I': 'I', 'SINE2/tRNA': 'tRNA',
                         'SINE1/7SL': '7SL', 'SINE3/5S': '5S', 'Unknown': 'Unknown'}

## Augmentation for each Repbase data
expandClassNum = {'Merlin': 20, 'Transib': 10, 'P': 10, 'Crypton': 10, 'Penelope': 5, 'R2': 20, 'RTE': 8, 'Jockey': 10, 'I': 10}

## Superfamilies classifiable by the ClassifyTE tool
ClassifyTE_class = {'Tc1-Mariner': '2.1.1.1', 'hAT': '2.1.1.2', 'Mutator': '2.1.1.3',
                'Merlin': '2.1.1.4', 'Transib': '2.1.1.5', 'P': '2.1.1.6',
                'PiggyBac': '2.1.1.7', 'PIF-Harbinger': '2.1.1.8', 'CACTA': '2.1.1.9',
                'Copia': '1.1.1', 'Gypsy': '1.1.2', 'Bel-Pao': '1.1.3',
                'DIRS': '1.2', 'R2': '1.4.1', 'RTE': '1.4.2',
                'Jockey': '1.4.3', 'L1': '1.4.4', 'I': '1.4.5',
                'tRNA': '1.5.1', '7SL': '1.5.2', '5S': '1.5.3'}

## Mapping Repbase labels to DeepTE labels
DeepTE_class = {'DNA_MITE_Tc': 'Tc1-Mariner', 'DNA_MITE_Harbinger': 'PIF-Harbinger',
             'DNA_MITE_hAT': 'hAT', 'DNA_MITE_CACTA': 'CACTA', 'DNA_MITE_MuDR': 'Mutator',
             'DNA_nMITE_Tc': 'Tc1-Mariner', 'DNA_nMITE_Harbinger': 'PIF-Harbinger', 'DNA_nMITE_hAT': 'hAT',
             'DNA_nMITE_CACTA': 'CACTA', 'DNA_nMITE_MuDR': 'Mutator', 'LTR_Copia': 'Copia', 'LTR_Gypsy': 'Gypsy',
             'LTR_ERV': 'Retrovirus', 'LTR_BEL': 'Bel-Pao', 'DIRS_DIRS': 'DIRS', 'RC_Helitron': 'Helitron'}
inverted_DeepTE_class = {'Tc1-Mariner': 'DNA_nMITE_Tc', 'PIF-Harbinger': 'DNA_nMITE_Harbinger', 'hAT': 'DNA_nMITE_hAT',
                         'CACTA': 'DNA_nMITE_CACTA', 'Mutator': 'DNA_nMITE_MuDR', 'Copia': 'LTR_Copia',
                         'Gypsy': 'LTR_Gypsy', 'Retrovirus': 'LTR_ERV', 'Bel-Pao': 'LTR_BEL', 'DIRS': 'DIRS_DIRS',
                         'Helitron': 'RC_Helitron'}