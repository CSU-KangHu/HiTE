import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader, SubsetRandomSampler
import pandas as pd
import numpy as np
torch.manual_seed(2024)
import torch.nn.functional as F
from mymodel import AttentionLSTMno
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
    return kmer_features

def build_counts_feature(char_array):
    counts = np.array([
    np.sum(char_array == 'A', axis=0),
    np.sum(char_array == 'C', axis=0),
    np.sum(char_array == 'G', axis=0),
    np.sum(char_array == 'T', axis=0),
    np.sum(char_array == '-', axis=0)
    ])
    return counts

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

def extract_features(file_path):
    ltr_names = []
    features = []
    for file in tqdm(file_path):
        ltr_names.append(file.split('.')[0].split('/')[-1])
        with open(file, 'r') as fr:
            lines = fr.readlines()
            char_array = np.array([list(s.replace('\t', '').strip()) for s in lines])
            counts_features = build_counts_feature(char_array)
            kmer_features = build_kmer_feature(char_array, [3,4,5])
            chunk_features = np.vstack((counts_features, kmer_features))
            features.extend(chunk_features)
    return features, ltr_names

def std_array(array):
    scaler = StandardScaler()
    return scaler.fit_transform(array)
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

def align2pileup(file_path, threads=20):
    # 进程数
    num_processes = threads  
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
    for result in results:
        chunk_feature, chunk_names = result
        features.extend(chunk_feature)
        ltr_names.extend(chunk_names)
    # print(features)
    features = np.array(features)
    std_features = std_array(features)
    return std_features, ltr_names

def main():
    # 从命令行接收参数
    parser = argparse.ArgumentParser(description="接收data_dir参数")
    parser.add_argument("--data_dir", type=str, help="数据目录路径")
    parser.add_argument("--out_dir", type=str, help="输出路径")    
    parser.add_argument("--model_path", type=str, help="模型路径")
    parser.add_argument("--threads", type=str, help="模型路径")
    args = parser.parse_args()
    outpath = args.out_dir
    file_path = args.data_dir
    model_path = args.model_path
    threads = int(args.threads)
    features, ltr_names = align2pileup(file_path, threads)
    
    result = []
    dim = 8
    for i in range(0, len(features), dim):
        chunk = features[i:i+dim]
        # combined = np.hstack(chunk)
        result.append(chunk)
    features = np.array(result) #（([100, 4, 100])
    # print(features[0])
    features = torch.tensor(features, dtype=torch.float32)
    # features = features.permute(0, 2, 1) # ([100, 100, 4])
    print(features.shape)
    dataset = torch.utils.data.TensorDataset(features)
    model = AttentionLSTMno()
    device = torch.device(f"cuda:{0}" if torch.cuda.is_available() else "cpu")
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()
    batch_size = 32
    val_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=threads)
    preds = []
    with torch.no_grad():
        for batch_x in tqdm(val_loader):
            # batch_x = batch_x.to(device)
            outputs = model(batch_x[0])
            _, predicted = torch.max(outputs, 1)
            preds.extend(predicted.cpu().numpy())
    with open(outpath + '/is_LTR.txt', 'w') as fw:
        for i in range(len(ltr_names)):
            fw.write(ltr_names[i] + '\t' + str(preds[i]) + '\n')

if __name__ == "__main__":
    main()