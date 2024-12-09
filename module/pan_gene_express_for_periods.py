#!/usr/bin/env python
import argparse
import os
import sys
import json
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, quantitative_gene


if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE detect de genes.")
    parser.add_argument("--genome_info", type=str, help="genome info json.")
    parser.add_argument("--RNA_dir", type=str, help="RNA sequence data directory.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    genome_info = args.genome_info
    RNA_dir = args.RNA_dir
    threads = args.threads

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/detect_de_genes.log', level='debug')

    # # Load the metadata
    # with open(genome_info, 'r') as f:
    #     genome_info_list = json.load(f)


    # # Step 7.2: 基因定量
    # log.logger.info("Start gene quantification using featureCounts...")
    # gene_express_dir = os.path.join(output_dir, 'gene_quantities')
    # os.makedirs(gene_express_dir, exist_ok=True)
    # RNA_tool_dir = os.path.join(project_dir, 'RNA_seq')
    # gene_express_table = quantitative_gene(genome_info_list, RNA_tool_dir, gene_express_dir, output_dir, threads, log)
    # log.logger.info(f"Gene quantification completed. Results saved to {gene_express_table}.")


    #gene_express_table = '/public/home/hpc194701009/ath_pan_genome/pan_genome/watermelon/panLTR_8_watermelon/G42_output/pan_gene_express_for_periods/gene_express.table'
    gene_express_table = '/home/hukang/test/HiTE/demo/gene_express.table'
    df = pd.read_csv(gene_express_table, sep="\t")
    df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: x.split(",")[-1])
    # print(df)

    samples = df.columns[1:]  # 除第一列之外的所有列
    x = np.arange(1, len(samples) + 1).reshape(-1, 1)  # 样本编号，作为自变量

    # 创建存储结果的列表
    results = []
    for _, row in df.iterrows():
        gene_name = row["gene_id"]
        y = row[samples].values.astype(float)
        mean = np.mean(y)
        std_dev = np.std(y)
        cv = std_dev / mean if mean != 0 else 0  # 避免除以0
        significant = cv > 0.5  # 自定义阈值
        results.append({"Gene_name": gene_name, "StdDev": std_dev, "CV": cv, "SignificantChange": significant})
    results_df = pd.DataFrame(results)
    print(results_df)
    results_df.to_csv("result.tsv", sep="\t", index=False)

    # 将差异基因和变化基因进行交集
    de_gene_table = '/home/hukang/test/HiTE/demo/DE_genes_from_TEs.tsv'
    de_gene_df = pd.read_csv(de_gene_table, sep="\t")
    print(de_gene_df)
    de_gene_df["Gene_name"] = de_gene_df["Gene_name"].str.replace('"', '')

    merged_df = pd.merge(results_df, de_gene_df, on="Gene_name", how="left")
    merged_df.to_csv("output.tsv", sep="\t", index=False)