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
random.seed(2024)

class FeatureExtractor:
    def __init__(self, input_dir, num_processes=10, logger=None):
        if logger is None:
            self.logger = self.make_logger()
        else:
            self.logger = logger
        self.input_dir = input_dir
        self.num_processes = num_processes
        self.target_files = self.get_files()
        self.Ks = [3, 4, 5]

    def cleanup(self):
        """
        clear all resources and manually call garbage collection.
        """
        self.logger.debug("Cleaning up resources...")

        # manually call garbage collection
        gc.collect()

        self.logger.debug("Resources cleaned up successfully.")

    def make_logger(self):
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
        process a single file and extract features
        """
        file_path = args
        try:
            with open(file_path, 'r') as fr:
                lines = fr.readlines()
                
            if len(lines) < 5:
                return None

            file_name = os.path.basename(file_path)
            ltr_name = file_name.split('.')[0].split('/')[-1]
            name = ltr_name
            
            lines = [re.sub(r'[^ATCG]', '-', line.replace('\t', '').strip()) for line in lines]
            char_array = np.array([list(s) for s in lines])
            char_array[char_array == 'N'] = '-'
            
            feature_row = 100
            padding_size = feature_row - char_array.shape[0]
            if padding_size > 0:
                pad_array = np.full((padding_size, char_array.shape[1]), '-', dtype=char_array.dtype)
                char_array = np.concatenate((char_array, pad_array), axis=0)
            elif padding_size < 0:
                selected_rows = random.sample(range(char_array.shape[0]), feature_row)
                char_array = char_array[selected_rows, :]
                
            char_array_left = char_array[:, :100]
            char_array_right = char_array[:, 100:]
            
            # extract features
            block_img_left = self.get_block_img(char_array_left)
            block_img_right = self.get_block_img(char_array_right)
            support_img_left = self.get_support_img(char_array_left)
            support_img_right = self.get_support_img(char_array_right)
            base_img_left = self.get_base_img(char_array_left)
            base_img_right = self.get_base_img(char_array_right)
            
            freq_left = self.get_freq_feature(char_array_left, self.Ks)
            freq_right = self.get_freq_feature(char_array_right, self.Ks)
            
            # merge features
            freq_feature = np.stack((freq_left, freq_right), axis=0)
            block_img = np.concatenate((block_img_left, block_img_right), axis=1)
            support_img = np.concatenate((support_img_left, support_img_right), axis=1)
            base_img = np.concatenate((base_img_left, base_img_right), axis=1)
            features = np.stack((block_img, support_img, base_img), axis=2)
            
            return {
                'features': features,
                'freqs': freq_feature,
                'name': name
            }
        except Exception as e:
            self.logger.error(f"Failed to process file {file_path}: {str(e)}")
            return None

    def generate_kmer_dic(self, repeat_num):
        kmer_dic = {}
        bases = ['A', 'G', 'C', 'T', '-']
        kmer_list = list(itertools.product(bases, repeat=int(repeat_num)))
        reduce_mer = '-' * repeat_num
        for eachitem in kmer_list:
            each_kmer = ''.join(eachitem)
            if each_kmer == reduce_mer:
                continue
            kmer_dic[each_kmer] = 0
        return kmer_dic

    def get_freq_feature(self, char_array, Ks):
        freq_features = []
        for k in Ks:
            kmer_dic = self.generate_kmer_dic(k)
            kmer_list_len = char_array.shape[1] - k + 1
            concatenated_strings = np.array([''.join(row[i:i+k]) for row in char_array for i in range(kmer_list_len)])
            concatenated_strings = concatenated_strings.reshape(-1, kmer_list_len)
            mers = []
            for i in range(kmer_list_len):
                mers += concatenated_strings[:, i].tolist()
            reduce_mer = '-' * k
            for mer in mers:
                if mer == reduce_mer:
                    continue
                kmer_dic[mer] += 1
            num_list = []
            for eachkmer in kmer_dic:
                num_list.append(kmer_dic[eachkmer])
            freq_features += num_list
        freq_features = np.array(freq_features)
        return freq_features

    def get_block_img(self, char_array):
        block_array = np.full(char_array.shape, 100, dtype=int)
        block_array[char_array == '-'] = 0
        return block_array

    def get_support_img(self, char_array):
        counts_array = np.full(char_array.shape, 100, dtype=int)
        letters = set('AGTC')
        for col in range(char_array.shape[1]):
            total_count = sum(np.array([np.sum(char_array[:, col] == letter) for letter in letters]))
            for letter in letters:
                if total_count == 0:
                    counts_array[char_array[:, col] == '-', col] = 0
                    break
                else:
                    count = np.sum(char_array[:, col] == letter)
                    counts_array[char_array[:, col] == letter, col] = int((count / total_count) * 100)
        counts_array[char_array == '-'] = 0
        return counts_array

    def get_base_img(self, char_array):
        base_dic = {'G': 100, 'C': 75, 'A': 50, 'T': 25, '-': 0}
        base_array = np.array([[base_dic[base] for base in row] for row in char_array])
        return base_array

    def get_files(self):
        start_time = time.time()
        self.logger.info("reading files...")
        pattern = os.path.join(self.input_dir, '*.matrix')
        target_files = glob.glob(pattern)
        end_time = time.time()
        self.logger.info(f"files read, time taken {end_time - start_time:.2f} seconds")

        return target_files

    def extract_main(self):

        self.logger.debug(f'number of samples, {len(self.target_files)}')

        # process files
        tasks = [f for f in self.target_files]
        self.logger.debug(f'number of tasks, {len(tasks)}')
        
        # use multiprocessing to process files
        results = []
        with Pool(processes=self.num_processes) as pool:
            for result in tqdm(pool.imap_unordered(self.process_file_chunk, tasks), 
                             total=len(tasks), 
                             desc="Processing files"):
                if result is not None:
                    results.append(result)

        # collect features
        features_list = []
        freqs_list = []
        seq_names = []
        
        for result in results:
            features_list.append(torch.from_numpy(result['features']))
            freqs_list.append(torch.from_numpy(result['freqs']))
            seq_names.append(result['name'])

        # merge features
        img_features = torch.stack(features_list)
        freq_features = torch.stack(freqs_list)
        
        # adjust dimensions
        img_features = img_features.permute(0, 3, 1, 2)
        freq_features = torch.unsqueeze(freq_features, dim=2)
        
        return img_features, freq_features, seq_names


def main():
    parser = argparse.ArgumentParser(description="Hybrid feature extraction program (Deep Learning)")
    parser.add_argument("--matrix_dir", type=str, required=True, help="LTR alignment matrix directory")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads")
    parser.add_argument("--feature_output_dir", type=str, required=True, help="Feature output directory")
    
    args = parser.parse_args()

    try:
        # create feature extractor instance
        feature_extractor = FeatureExtractor(args.matrix_dir, args.threads)

        if not os.path.exists(args.feature_output_dir):
            os.makedirs(args.feature_output_dir, exist_ok=True)
            feature_extractor.logger.info(f"Feature output folder -- {args.feature_output_dir} created successfully")
        
        # extract features and save results
        img_features, freq_features, seq_names = feature_extractor.extract_main()
        torch.save(img_features, os.path.join(args.feature_output_dir, 'img_features.pt'))
        torch.save(freq_features, os.path.join(args.feature_output_dir, 'freq_features.pt'))
        pd.DataFrame(seq_names).to_csv(os.path.join(args.feature_output_dir, 'seq_names.txt'), index=False, header=False)

        feature_extractor.cleanup()
    except Exception as e:
        logging.error(f"Feature extraction failed: {str(e)}")
        raise
