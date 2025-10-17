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
is_debug = 0 # Is debug mode


# 2. Program and model parameters
kmer_sizes = [1, 3, 5]
threads = int(cpu_count())  # Use multi-threading to load data
## CNN model parameters
cnn_num_convs = 5 # Number of CNN convolutional layers
# cnn_filters_array = [32, 64, 64] # Number of filters per convolutional layer in CNN
# cnn_kernel_sizes_array = [(3, 3), (3, 3), (3, 3)] # Kernel size for each convolutional layer in CNN; for 2D convolutional layers, set as [(3, 3), ...]
cnn_filters_array = [32, 32, 32, 32, 32] # Number of filters per convolutional layer in CNN
cnn_kernel_sizes_array = [3, 3, 3, 3, 3] # Kernel size for each convolutional layer in CNN; for 2D convolutional layers, set as [(3, 3), ...]
cnn_dropout = 0.5 # CNN dropout threshold
## Training parameters
batch_size = 32 # Batch size for training
epochs = 25 #  Number of epochs for training
use_checkpoint = 0  # Whether to use checkpoint training; set to 1 to resume training from the parameters of the previous failed training session, avoiding training from scratch

lstm_num_layers = 1
lstm_hidden_size = 128
lstm_dropout = 0.5

################################################### The following parameters do not need modification ######################################################################
version_num = '1.0.1'
work_dir = project_dir + '/work' # temp work directory

X_feature_len = 0
# for kmer_size in kmer_sizes:
#     X_feature_len += 100 - kmer_size + 1
for kmer_size in kmer_sizes:
    X_feature_len += pow(5, kmer_size)
class_num = 2

input_shape = (X_feature_len, 1)

# Obtain CNN input dimensions
image_shape = (100, 200, 3)

# Obtain LSTM input dimensions
lstm_seq_len = 100
lstm_features_num = 5