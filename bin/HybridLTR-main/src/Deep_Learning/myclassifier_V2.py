import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader, SubsetRandomSampler
import pandas as pd
import numpy as np
torch.manual_seed(2024)
import torch.nn.functional as F
from mymodel_V2 import LSTMCat
import os
from tqdm import tqdm
import sys
import argparse
import math
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances
import itertools
import re
import multiprocessing


def chunk_data(pos_data, chunk_size1):
    for i in range(0, len(pos_data), chunk_size1):
        pos_chunk = pos_data[i:i + chunk_size1] if i is not None else []
        yield pos_chunk

def build_kmer_feature(char_array, Ks):
    kmer_features = []
    for k in Ks:
        kmer_list_len = char_array.shape[1] - k + 1
        concatenated_strings = np.array([''.join(row[i:i+k]) for row in char_array for i in range(kmer_list_len)])
        concatenated_strings = concatenated_strings.reshape(-1, kmer_list_len)
        counts_mer = [len(set(concatenated_strings[:,i])) for i in range(kmer_list_len)]  
        counts_mer = counts_mer + [0] * (k-1)
        kmer_features.append(counts_mer)
    kmer_features = np.array(kmer_features)
    column_sums = kmer_features.sum(axis=0)
    zero_columns = column_sums == 0
    kmer_ratios = np.zeros_like(kmer_features)
    for i in range(kmer_features.shape[1]):
        if not zero_columns[i]:  # 如果列的总和不为0
            kmer_ratios[:, i] = kmer_features[:, i] / column_sums[i]

    return kmer_ratios

def build_counts_feature(char_array):
    counts = np.array([
    np.sum(char_array == 'A', axis=0),
    np.sum(char_array == 'C', axis=0),
    np.sum(char_array == 'G', axis=0),
    np.sum(char_array == 'T', axis=0),
    np.sum(char_array == '-', axis=0)
    ])
    column_sums = counts.sum(axis=0)
    counts_ratios = counts / column_sums
    return counts_ratios

def generate_kmer_dic(repeat_num):
    ##initiate a dic to store the kmer dic
    ##kmer_dic = {'ATC':0,'TTC':1,...}
    kmer_dic = {}
    bases = ['A','G','C','T', '-']
    kmer_list = list(itertools.product(bases, repeat=int(repeat_num)))
    for eachitem in kmer_list:
        each_kmer = ''.join(eachitem)
        kmer_dic[each_kmer] = 0
    return (kmer_dic)

def build_freq_feature(char_array, Ks):
    freq_features = []
    for k in Ks:
        kmer_dic = generate_kmer_dic(k)
        kmer_list_len = char_array.shape[1] - k + 1
        concatenated_strings = np.array([''.join(row[i:i+k]) for row in char_array for i in range(kmer_list_len)])
        concatenated_strings = concatenated_strings.reshape(-1, kmer_list_len)
        mers = []
        for i in range(kmer_list_len):
            mers += concatenated_strings[:,i].tolist()
        for mer in mers:
            kmer_dic[mer] += 1
        num_list = []  ##this dic stores num_dic = [0,1,1,0,3,4,5,8,2...]
        for eachkmer in kmer_dic:
            num_list.append(kmer_dic[eachkmer])
        freq_features += num_list
    freq_features = np.array(freq_features)
    total_sum = freq_features.sum()
    freq_ratio = freq_features / total_sum
    return freq_ratio

def extract_features(file_path):
    ltr_names = []
    features = []
    chunk_freq_features = []

    for file in tqdm(file_path):
        ltr_names.append(file.split('.')[0].split('/')[-1])
        with open(file, 'r') as fr:
            lines = fr.readlines()
            lines = [re.sub(r'[^ATCG]', '-', line.replace('\t', '').strip()) for line in lines]
            char_array = np.array([list(s) for s in lines])
            counts_features = build_counts_feature(char_array)
            kmer_features = build_kmer_feature(char_array, [3,4,5])
            chunk_features = np.vstack((counts_features, kmer_features))
            features.extend(chunk_features)
            left_half = np.array([row[:len(row)//2] for row in char_array])
            freq_features_left = build_freq_feature(left_half, [2,3,4])
            right_half = np.array([row[len(row)//2:] for row in char_array])
            freq_features_right = build_freq_feature(right_half, [2,3,4])
            chunk_freq_features.append(freq_features_left)
            chunk_freq_features.append(freq_features_right)
    
    return features, chunk_freq_features, ltr_names

def std_array(array):
    scaler = StandardScaler()
    return scaler.fit_transform(array)

def transform_freqarray(freq_features):
    leftarray = [freq_features[i] for i in range(0, len(freq_features), 2)]
    rightarray = [freq_features[i] for i in range(1, len(freq_features), 2)]
    cat_freq_features = np.stack((np.array(leftarray), np.array(rightarray)), axis=1)
    return cat_freq_features

# 读取文件
def get_files(file_path):
    target_files = []
    for root, dirs, files in os.walk(file_path):
        for file in files:
            # 检查文件是否以'.matrix'结尾
            if file.endswith('.matrix'):
                # 构建完整路径
                file_name = os.path.join(root, file)
                # 存储路径
                target_files.append(file_name)
    return target_files

def align2pileup(file_path):
    # 进程数
    num_processes = 20  
    target_files = get_files(file_path)
    print(len(target_files))
    chunk_size = len(target_files) // num_processes 
    data_chunks = list(chunk_data(target_files, chunk_size))
    # 创建进程池
    # pool = multiprocessing.Pool(processes=num_processes)
    # results = pool.map(extract_features, data_chunks)
    with multiprocessing.Pool(processes=num_processes) as pool:
        # 使用 map 方法将 process_array 函数应用于每个数组
        results = pool.map(extract_features,  data_chunks)
    features = []
    ltr_names = []
    freq_features = []
    for result in results:
        chunk_feature, chunk_freq_feature, chunk_names = result
        features.extend(chunk_feature)
        freq_features.extend(chunk_freq_feature)
        ltr_names.extend(chunk_names)
    # print(features)
    features = np.array(features)
    std_features = std_array(features)
    freq_features = np.array(freq_features)
    std_freq_features = freq_features
    std_freq_features = std_array(freq_features) # 暂时不标准化
    return std_features, std_freq_features, ltr_names

def main():
    # 从命令行接收参数
    parser = argparse.ArgumentParser(description="接收data_dir参数")
    parser.add_argument("--data_dir", type=str, help="数据目录路径")
    parser.add_argument("--out_dir", type=str, help="输出路径")    
    parser.add_argument("--model_path", type=str, help="模型路径")
    args = parser.parse_args()
    outpath = args.out_dir
    file_path = args.data_dir
    model_path = args.model_path

    features, freq_features, ltr_names = align2pileup(file_path)
    
    result = []
    dim = 8
    for i in range(0, len(features), dim):
        chunk = features[i:i+dim]
        # combined = np.hstack(chunk)
        result.append(chunk)
    features = np.array(result) #（([100, 4, 100])
    # print(features[0])
    features = torch.tensor(features, dtype=torch.float32)
    freq_features = transform_freqarray(freq_features)
    # print(freq_features[0])
    freq_features = torch.tensor(freq_features, dtype=torch.float32)
    # features = features.permute(0, 2, 1) # ([100, 100, 4])
    print(features.shape)
    print(freq_features.shape)
    dataset = torch.utils.data.TensorDataset(features, freq_features)
    model = LSTMCat()
    device = torch.device(f"cuda:{0}" if torch.cuda.is_available() else "cpu")
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()
    batch_size = 32
    val_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    preds = []
    with torch.no_grad():
        for batch_x, bathc_y in tqdm(val_loader):
            outputs = model(batch_x, bathc_y)
            _, predicted = torch.max(outputs, 1)
            preds.extend(predicted.cpu().numpy())
    with open(outpath + '/is_LTR.txt', 'w') as fw:
        for i in range(len(ltr_names)):
            fw.write(ltr_names[i] + '\t' + str(preds[i]) + '\n')

if __name__ == "__main__":
    main()