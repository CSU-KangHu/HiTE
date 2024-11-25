import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader, SubsetRandomSampler
import pandas as pd
import numpy as np
torch.manual_seed(2024)
import torch.nn.functional as F
from mymodel import AttentionLSTMCNN, AttentionLSTMCNN3, AttentionLSTMno, AttentionLSTMCNN4
import os
from tqdm import tqdm
import sys
import argparse
import math
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances

def log_transform(x, base):
    x = int(x)
    if x == 0:
        return 0
    else:
       return math.log(x, base)

def transform_array(arr):
    transformed_arr = []
    for row in arr:
        transformed_row = [log_transform(value, 10) for value in row]
        transformed_arr.append(transformed_row)
    return transformed_arr


def hamming_distance(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

# 计算汉明距离矩阵
def hamming_distance_matrix(data):
    n = len(data)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            distance_matrix[i, j] = hamming_distance(data[i], data[j])
    return distance_matrix

def build_kmer_featureV2(char_array, Ks):
    kmer_features = []
    for k in Ks:
        kmer_list_len = char_array.shape[1] - k + 1
        concatenated_strings = np.array([''.join(row[i:i+k]) for row in char_array for i in range(kmer_list_len)])
        concatenated_strings = concatenated_strings.reshape(-1, kmer_list_len)
        counts_mer = []
        for i in range(kmer_list_len):
            mers = concatenated_strings[:,i]
            distance_matrix = hamming_distance_matrix(mers)
            dbscan = DBSCAN(metric="precomputed", eps=1, min_samples=1)
            labels = dbscan.fit_predict(distance_matrix)
            counts_mer.append(len(set(labels)))
        counts_mer = counts_mer + [0] * (k-1)
        kmer_features.append(counts_mer)
    return np.array(kmer_features)

def build_kmer_feature(char_array, Ks):
    kmer_features = []
    for k in Ks:
        kmer_list_len = char_array.shape[1] - k + 1
        concatenated_strings = np.array([''.join(row[i:i+k]) for row in char_array for i in range(kmer_list_len)])
        concatenated_strings = concatenated_strings.reshape(-1, kmer_list_len)
        counts_mer = [len(set(concatenated_strings[:,i])) for i in range(kmer_list_len)]  
        counts_mer = counts_mer + [0] * (k-1)
        # counts_mer_log = [np.log10(count) if count > 0 else 0 for count in counts_mer]
        # counts_mer_log = counts_mer
        # kmer_features.append(counts_mer_log)
        kmer_features.append(counts_mer)
    kmer_features = np.array(kmer_features)
    # scaler = MinMaxScaler()
    # kmer_features = scaler.fit_transform(kmer_features)
    return kmer_features

def build_counts_feature(char_array):
    counts = np.array([
    np.sum(char_array == 'A', axis=0),
    np.sum(char_array == 'C', axis=0),
    np.sum(char_array == 'G', axis=0),
    np.sum(char_array == 'T', axis=0),
    np.sum(char_array == '-', axis=0)
    ])
    # counts_log = np.log1p(counts)
    # counts_log = counts
    # counts_log = np.where(counts > 0, np.log(counts), counts)
    # scaler = MinMaxScaler()
    # normalized_array = scaler.fit_transform(counts)
    return counts


def extract_features(file_path):
    ltr_files = []
    ltr_names = []
    for root, dirs, files in os.walk(file_path):
        for file in files:
            # 检查文件是否以'.matrix'结尾
            if file.endswith('.matrix'):
                ltr_names.append(file.split('.')[0])
                # 构建完整路径
                file_path = os.path.join(root, file)
                # 存储路径
                ltr_files.append(file_path)
    features = []
    for file in tqdm(ltr_files):
        with open(file, 'r') as fr:
            lines = fr.readlines()
            char_array = np.array([list(s.replace('\t', '').strip()) for s in lines])
            counts_features = build_counts_feature(char_array)
            kmer_features = build_kmer_feature(char_array, [3,4,5])
            chunk_features = np.vstack((counts_features, kmer_features))
            features.extend(chunk_features)
            # temp_feature = [np.zeros(100) for _ in range(5)]
            # for line in lines:
            #     line = line.strip()
            #     for char_idx in range(len(line)):
            #         if line[char_idx] == 'A':
            #             temp_feature[0][char_idx] += 1
            #         elif line[char_idx] == 'G':
            #             temp_feature[1][char_idx] += 1
            #         elif line[char_idx] == 'C':
            #             temp_feature[2][char_idx] += 1
            #         elif line[char_idx] == 'T':
            #             temp_feature[3][char_idx] += 1
            #         elif line[char_idx] == '-':
            #             temp_feature[4][char_idx] += 1
            # temp_feature_ratio = [np.zeros(100) for _ in range(4)]
            # sums = [sum(column[:4]) for column in zip(*temp_feature)]
            # for i in range(100):
            #     for j in range(4):
            #         if sums[i] != 0:
            #             temp_feature_ratio[j][i] = temp_feature[j][i] / sums[i]
            #         else:
            #             temp_feature_ratio[j][i] = 0
            # features.extend(temp_feature_ratio + temp_feature)
            # features.extend(transform_array(temp_feature))
                        
    return features, ltr_names

def std_array(array):
    scaler = StandardScaler()
    return scaler.fit_transform(array)

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
    features, ltr_names = extract_features(file_path)
    features = std_array(features)
    result = []
    dim = 8
    for i in range(0, len(features), dim):
        chunk = features[i:i+dim]
        # combined = np.hstack(chunk)
        result.append(chunk)
    features = np.array(result) #（([100, 4, 100])
    features = torch.tensor(features, dtype=torch.float32)
    # features = features.permute(0, 2, 1) # ([100, 100, 4])
    print(features.shape)
    dataset = torch.utils.data.TensorDataset(features)
    model = AttentionLSTMno()
    device = torch.device(f"cuda:{0}" if torch.cuda.is_available() else "cpu")
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()
    batch_size = 32
    val_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    preds = []
    with torch.no_grad():
        for batch_x in tqdm(val_loader):
            outputs = model(batch_x[0])
            _, predicted = torch.max(outputs, 1)
            preds.extend(predicted.cpu().numpy())
    with open(outpath + '/is_LTR.txt', 'w') as fw:
        for i in range(len(ltr_names)):
            fw.write(ltr_names[i] + '\t' + str(preds[i]) + '\n')

if __name__ == "__main__":
    main()