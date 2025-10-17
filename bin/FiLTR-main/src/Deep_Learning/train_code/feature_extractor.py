import os
import re
import math
import glob
import random
import numpy as np
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm
import multiprocessing
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import itertools
from itertools import zip_longest
import torch
import torch.nn.functional as F
import logging
import time
import gc

# Set a seed for reproducibility
random.seed(2024)

class FeatureExtractor:
    """
    Extracts features from sequence alignment files (.matrix) for model training.
    It processes files in parallel, generates image-like and frequency-based features,
    and saves them as PyTorch tensors.
    """
    def __init__(self, input_dir, num_processes=10):
        self.logger = self.make_logger()
        self.input_dir = input_dir
        self.num_processes = num_processes
        self.pos_files, self.neg_files = self.get_files()
        self.Ks = [3, 4, 5] # k-mer sizes for frequency feature extraction

    def cleanup(self):
        """
        Cleans up all resources and manually triggers garbage collection.
        """
        print("Cleaning up resources...")
        # Manually trigger garbage collection
        gc.collect()
        print("Resources cleaned up successfully.")

    def make_logger(self):
        """Initializes a console logger."""
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        return logger

    def process_file_chunk(self, args):
        """
        Processes a single file to extract features.
        
        Args:
            args (tuple): A tuple containing the file path and its corresponding label.
            
        Returns:
            A dictionary containing the extracted features, frequency data, label, and a unique name.
        """
        file_path, label_value = args
        try:
            with open(file_path, 'r') as fr:
                lines = fr.readlines()
                
            # Skip files with insufficient data
            if len(lines) < 5:
                return None
                
            # Create a unique name for the sequence
            ltr_name = file_path.split('.')[0].split('/')[-1]
            name = f"{label_value}_{ltr_name}%{len(lines)}"
            
            # Clean and format the sequence data
            lines = [re.sub(r'[^ATCG]', '-', line.replace('\t', '').strip()) for line in lines]
            char_array = np.array([list(s) for s in lines])
            char_array[char_array == 'N'] = '-'
            
            # Pad or sample rows to a fixed height of 100
            feature_row = 100
            padding_size = feature_row - char_array.shape[0]
            if padding_size > 0:
                pad_array = np.full((padding_size, char_array.shape[1]), '-', dtype=char_array.dtype)
                char_array = np.concatenate((char_array, pad_array), axis=0)
            elif padding_size < 0:
                selected_rows = random.sample(range(char_array.shape[0]), feature_row)
                char_array = char_array[selected_rows, :]
                
            # Split the array into left and right halves
            char_array_left = char_array[:, :100]
            char_array_right = char_array[:, 100:]
            
            # Extract image-like features
            block_img_left = self.get_block_img(char_array_left)
            block_img_right = self.get_block_img(char_array_right)
            support_img_left = self.get_support_img(char_array_left)
            support_img_right = self.get_support_img(char_array_right)
            base_img_left = self.get_base_img(char_array_left)
            base_img_right = self.get_base_img(char_array_right)
            
            # Extract k-mer frequency features
            freq_left = self.get_freq_feature(char_array_left, self.Ks)
            freq_right = self.get_freq_feature(char_array_right, self.Ks)
            
            # Assemble the final features
            freq_feature = np.stack((freq_left, freq_right), axis=0)
            block_img = np.concatenate((block_img_left, block_img_right), axis=1)
            support_img = np.concatenate((support_img_left, support_img_right), axis=1)
            base_img = np.concatenate((base_img_left, base_img_right), axis=1)
            features = np.stack((block_img, support_img, base_img), axis=2)
            
            return {
                'features': features,
                'freqs': freq_feature,
                'label': label_value,
                'name': name
            }
        except Exception as e:
            self.logger.error(f"Error processing file {file_path}: {e}")
            return None

    def generate_kmer_dic(self, k):
        """Generates a dictionary for all possible k-mers with initial counts of 0."""
        kmer_dic = {}
        bases = ['A', 'G', 'C', 'T', '-']
        kmer_list = list(itertools.product(bases, repeat=int(k)))
        reduce_mer = '-' * k
        for item in kmer_list:
            each_kmer = ''.join(item)
            if each_kmer == reduce_mer:
                continue
            kmer_dic[each_kmer] = 0
        return kmer_dic

    def get_freq_feature(self, char_array, Ks):
        """Calculates k-mer frequencies for a given character array."""
        freq_features = []
        for k in Ks:
            kmer_dic = self.generate_kmer_dic(k)
            num_kmers_per_row = char_array.shape[1] - k + 1
            
            # Efficiently extract all k-mers
            all_kmers = [''.join(row[i:i+k]) for row in char_array for i in range(num_kmers_per_row)]
            
            # Count k-mer occurrences
            reduce_mer = '-' * k
            for mer in all_kmers:
                if mer != reduce_mer and mer in kmer_dic:
                    kmer_dic[mer] += 1
            
            freq_features.extend(list(kmer_dic.values()))
            
        return np.array(freq_features)

    def get_block_img(self, char_array):
        """Creates a binary mask where bases are 100 and gaps ('-') are 0."""
        block_array = np.full(char_array.shape, 100, dtype=int)
        block_array[char_array == '-'] = 0
        return block_array

    def get_support_img(self, char_array):
        """Creates a support map based on the frequency of each base in its column."""
        counts_array = np.full(char_array.shape, 100, dtype=int)
        letters = set('AGTC')
        for col in range(char_array.shape[1]):
            column_data = char_array[:, col]
            total_count = np.sum(np.isin(column_data, list(letters)))
            
            if total_count == 0:
                counts_array[column_data == '-', col] = 0
                continue
            
            for letter in letters:
                count = np.sum(column_data == letter)
                counts_array[column_data == letter, col] = int((count / total_count) * 100)
                
        counts_array[char_array == '-'] = 0
        return counts_array

    def get_base_img(self, char_array):
        """Creates a numeric representation of the alignment based on a fixed mapping."""
        base_dic = {'G': 100, 'C': 75, 'A': 50, 'T': 25, '-': 0}
        base_array = np.vectorize(base_dic.get)(char_array)
        return base_array

    def get_files(self):
        """Retrieves the list of positive and negative sample files from the input directory."""
        start_time = time.time()
        self.logger.info("Reading files...")
        neg_path = os.path.join(self.input_dir, "negative")
        pos_path = os.path.join(self.input_dir, "positive")
        neg_files = glob.glob(os.path.join(neg_path, '*.matrix'))
        pos_files = glob.glob(os.path.join(pos_path, '*.matrix'))
        self.logger.debug(f'Positive files found: {len(pos_files)}')
        self.logger.debug(f'Negative files found: {len(neg_files)}')
        end_time = time.time()
        self.logger.info(f"File reading completed in {end_time - start_time:.2f} seconds")
        return pos_files, neg_files

    def main(self):
        """Main execution function to process all files and return feature tensors."""
        print('Positive samples:', len(self.pos_files))
        print('Negative samples:', len(self.neg_files))

        # Prepare a list of tasks for the process pool
        tasks = [(f, 1) for f in self.pos_files] + [(f, 0) for f in self.neg_files]
        
        # Use a process pool to handle tasks in parallel
        results = []
        with Pool(processes=self.num_processes) as pool:
            # Use tqdm for a progress bar
            for result in tqdm(pool.imap_unordered(self.process_file_chunk, tasks), 
                             total=len(tasks), 
                             desc="Processing files"):
                if result is not None:
                    results.append(result)

        # Collect and organize the results
        features_list = [torch.from_numpy(r['features']) for r in results]
        freqs_list = [torch.from_numpy(r['freqs']) for r in results]
        labels = [r['label'] for r in results]
        seq_names = [r['name'] for r in results]
        
        # Stack features into tensors
        img_features = torch.stack(features_list)
        freq_features = torch.stack(freqs_list)
        
        # Adjust tensor dimensions for the model input
        # (N, H, W, C) -> (N, C, H, W)
        img_features = img_features.permute(0, 3, 1, 2)
        # Add a channel dimension for frequency features
        freq_features = torch.unsqueeze(freq_features, dim=2)
        
        # Convert labels to a tensor
        labels = torch.tensor(labels, dtype=torch.long)

        return img_features, freq_features, labels, seq_names


if __name__ == '__main__':
    # Define directories and number of processes
    train_input_dir = "your/train/data/path/"
    test_input_dir = "your/test/data/path/"
    num_processes = 20

    # Process training data
    print("--- Processing Training Data ---")
    feature_extractor_train = FeatureExtractor(train_input_dir, num_processes)
    train_img_features, train_freq_features, train_labels, train_seq_names = feature_extractor_train.main()
    
    # Save training tensors and names
    torch.save(train_img_features, 'train_img_features.pt')
    torch.save(train_freq_features, 'train_freq_features.pt')
    torch.save(train_labels, 'train_labels.pt')
    with open('train_seq_names.txt', 'w') as fw:
        for name in train_seq_names:
            fw.write(name + '\n')
            
    print('train_img_features shape:', train_img_features.shape)
    print('train_freq_features shape:', train_freq_features.shape)
    print('train_labels shape:', train_labels.shape)
    feature_extractor_train.cleanup()

    # Process testing data
    print("\n--- Processing Test Data ---")
    feature_extractor_test = FeatureExtractor(test_input_dir, num_processes)
    test_img_features, test_freq_features, test_labels, test_seq_names = feature_extractor_test.main()
    
    # Save testing tensors and names
    torch.save(test_img_features, 'test_img_features.pt')
    torch.save(test_freq_features, 'test_freq_features.pt')
    torch.save(test_labels, 'test_labels.pt')
    with open('test_seq_names.txt', 'w') as fw:
        for name in test_seq_names:
            fw.write(name + '\n')
            
    print('test_img_features shape:', test_img_features.shape)
    print('test_freq_features shape:', test_freq_features.shape)
    print('test_labels shape:', test_labels.shape)
    feature_extractor_test.cleanup()