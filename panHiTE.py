import argparse
import json
import os
import sys
import time
from multiprocessing import cpu_count
import subprocess

current_folder = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.join(current_folder, ".")

from module.Util import read_fasta, Logger, \
    store_fasta, lib_add_prefix, find_gene_relation_tes, deredundant_for_LTR_v5, \
    generate_bam_for_RNA_seq, quantitative_gene, summary_TEs, ReassignInconsistentLabels, file_exist, \
    convertGeneAnnotation2GTF


def main_pipeline():
    tool_name = 'panHiTE_serial'
    version_num = '0.0.1'
    default_threads = int(cpu_count())
    default_miu = str(1.3e-8)
    default_recover = 0
    default_te_type = 'all'
    default_skip_analyze = 0
    default_debug = 0

    # 1.parse args
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--pan_genomes_dir', required=True, metavar='pan_genomes_dir', help='A directory containing the pangenome.')
    parser.add_argument('--genome_list', required=True, metavar='genome_list', help='A text file with genome and gene names. Each line represents a pair of genome and gene names, separated by a tab (\t). The genome name is mandatory, while the gene name is optional. If a gene name is provided, the genes_dir parameter must also be specified.')
    parser.add_argument('--genes_dir', required=False, metavar='genes_dir', help='A directory containing the sorted gene annotation files, gff or gtf format.')
    parser.add_argument('--RNA_dir', required=False, metavar='RNA_dir', help='A directory containing the RNA-seq files.')

    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.')

    parser.add_argument('--thread', metavar='thread_num', help='Input thread num, default = [ ' + str(default_threads) + ' ]')
    parser.add_argument('--miu', metavar='miu', help='The neutral mutation rate (per bp per ya), default = [ ' + str(default_miu) + ' ]')
    parser.add_argument('--te_type', metavar='te_type', help='Retrieve specific type of TE output [ltr|tir|helitron|non-ltr|all]. default = [ ' + str(default_te_type) + ' ]')
    parser.add_argument('--skip_analyze', metavar='skip_analyze', help='Whether to skip analyze, only generate panTE library. default = [ ' + str(default_skip_analyze) + ' ]')
    parser.add_argument('--debug', metavar='is_debug', help='Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ ' + str(default_debug) + ' ]')
    parser.add_argument('--recover', metavar='is_recover', help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ ' + str(default_recover) + ' ]')

    args = parser.parse_args()

    pan_genomes_dir = args.pan_genomes_dir
    genome_list_path = args.genome_list
    genes_dir = args.genes_dir
    RNA_dir = args.RNA_dir
    output_dir = args.out_dir
    threads = args.thread
    miu = args.miu
    te_type = args.te_type
    skip_analyze = args.skip_analyze
    recover = args.recover
    debug = args.debug

    output_dir = os.path.abspath(output_dir + '/')
    if os.path.exists(output_dir) and not recover:
        os.system('rm -rf ' + output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    if miu is None:
        miu = default_miu
    else:
        miu = str(miu)

    if te_type is None:
        te_type = default_te_type
    else:
        te_type = te_type

    if debug is None:
        debug = default_debug
    else:
        debug = int(debug)

    if skip_analyze is None:
        skip_analyze = default_skip_analyze
    else:
        skip_analyze = int(skip_analyze)

    if recover is None:
        recover = default_recover
    else:
        recover = int(recover)

    all_te_types = ['ltr', 'tir', 'helitron', 'non-ltr', 'all']
    if te_type not in all_te_types:
        log.logger.error('Specified an invalid TE type: ' + te_type + '. Please choose from ' + str(all_te_types))
        sys.exit(-1)

    # 1. 预处理步骤
    output_metadata = os.path.join(output_dir, 'genome_metadata.json')
    resut_file = output_metadata
    if not recover or not file_exist(resut_file):
        starttime = time.time()
        pan_preprocess_cmd = 'pan_preprocess_genomes.py ' + genome_list_path + ' ' + genes_dir \
                             + ' ' + RNA_dir + ' ' + pan_genomes_dir + ' ' + output_dir
        log.logger.info(pan_preprocess_cmd)
        os.system(pan_preprocess_cmd)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step1: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')

    # Load the metadata
    with open(output_metadata, 'r') as f:
        genome_data = json.load(f)
    total_genome = genome_data["total_genome"]
    genome_info_list = genome_data["genome_info"]
    gene_annotation_list = genome_data["gene_annotations"]
    is_RNA_analyze = genome_data["is_RNA_analyze"]


    # 2. 遍历每个基因组，并调用HiTE进行TE识别
    intact_ltr_paths = []
    pan_terminal_tmp_lib = os.path.join(output_dir, 'pan_terminal.tmp.fa')
    pan_internal_tmp_lib = os.path.join(output_dir, 'pan_internal.tmp.fa')
    panTE_lib = os.path.join(output_dir, 'panTE.fa')

    resut_file = panTE_lib
    if not recover or not file_exist(resut_file):
        starttime = time.time()
        for genome_info in genome_info_list:
            genome_name = genome_info["genome_name"]
            raw_name = genome_info["raw_name"]
            reference = genome_info["reference"]

            TE_gff = genome_info["TE_gff"]
            gene_gtf = genome_info["gene_gtf"]
            RNA_seq_dict = genome_info["RNA_seq"]

            HiTE_output_dir = output_dir + '/HiTE_' + str(raw_name)
            single_result_path = os.path.join(output_dir, f"{genome_name}_hite_result.json")
            pan_run_hite_single_cmd = 'pan_run_hite_single.py ' + genome_name + ' ' + reference \
                                      + ' ' + HiTE_output_dir + ' ' + str(threads) + ' ' + te_type \
                                      + ' ' + str(miu) +' ' + str(debug) + ' ' + str(recover)
            log.logger.info(pan_run_hite_single_cmd)
            os.system(pan_run_hite_single_cmd)

            # Load the metadata
            with open(single_result_path, 'r') as f:
                single_result = json.load(f)

            confident_ltr_terminal = single_result["confident_ltr_terminal"]
            confident_ltr_internal = single_result["confident_ltr_internal"]
            confident_TE = single_result["confident_TE"]
            ltr_intact_list = single_result["ltr_intact_list"]
            intact_ltr_paths.append(ltr_intact_list)
            confident_helitron = single_result["confident_helitron"]
            confident_non_ltr = single_result["confident_non_ltr"]
            confident_other = single_result["confident_other"]
            confident_tir = single_result["confident_tir"]

            # 2.3 将每个基因组生成的library合并
            os.system('cat ' + confident_ltr_terminal + ' >> ' + pan_terminal_tmp_lib)
            os.system('cat ' + confident_ltr_internal + ' >> ' + pan_internal_tmp_lib)
            os.system('cat ' + confident_TE + ' >> ' + pan_terminal_tmp_lib)

        pan_remove_redundancy_cmd = 'pan_remove_redundancy.py ' + pan_terminal_tmp_lib + ' ' + pan_internal_tmp_lib \
                                  + ' ' + output_dir + ' ' + str(threads) + ' ' + panTE_lib
        log.logger.info(pan_remove_redundancy_cmd)
        os.system(pan_remove_redundancy_cmd)

        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step1: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')
        for genome_info in genome_info_list:
            raw_name = genome_info["raw_name"]
            HiTE_output_dir = output_dir + '/HiTE_' + str(raw_name)
            ltr_intact_list = HiTE_output_dir + '/intact_LTR.list'
            intact_ltr_paths.append(ltr_intact_list)

    batch_files = []
    if not skip_analyze:
        # 4. 利用 panTE 注释每个基因组
        for genome_info in genome_info_list:
            genome_name = genome_info["genome_name"]
            raw_name = genome_info["raw_name"]
            reference = genome_info["reference"]

            TE_gff = genome_info["TE_gff"]
            gene_gtf = genome_info["gene_gtf"]
            RNA_seq_dict = genome_info["RNA_seq"]

            HiTE_output_dir = output_dir + '/HiTE_' + str(raw_name)
            single_result_path = os.path.join(output_dir, f"{genome_name}_hite_result.json")
            pan_run_hite_single_cmd = 'pan_run_hite_single.py ' + genome_name + ' ' + reference \
                                      + ' ' + HiTE_output_dir + ' ' + str(threads) + ' ' + te_type \
                                      + ' ' + str(miu) + ' ' + str(debug) + ' ' + str(recover)
            log.logger.info(pan_run_hite_single_cmd)
            os.system(pan_run_hite_single_cmd)



        for genome_name, reference, TE_gff, gene_gtf, RNA_seq_dict in genome_paths:
            full_length_TE_gff = output_dir + '/' + genome_name + '.full_length.gff'
            batch_files.append((genome_name, reference, TE_gff, full_length_TE_gff, gene_gtf, RNA_seq_dict))
            resut_file = TE_gff
            if not recover or not os.path.exists(resut_file):
                RepeatMasker_command = 'cd ' + output_dir + ' && RepeatMasker -e ncbi -no_is -norna -nolow -pa ' + str(
                    threads) + ' -gff -lib ' + panTE_lib + ' -cutoff 225 ' + reference
                os.system(RepeatMasker_command)

                mv_file_command = 'mv ' + reference + '.out ' + output_dir + '/' + genome_name + '.out && mv ' \
                                  + reference + '.tbl ' + output_dir + '/' + genome_name + '.tbl && mv ' \
                                  + reference + '.out.gff ' + output_dir + '/' + genome_name + '.gff'
                os.system(mv_file_command)
            else:
                log.logger.info(resut_file + ' exists, skip...')


        # 5. 根据注释好的 annotation 文件，进行常见TE分析
        log.logger.info('Start analysing using TE annotation files...')
        summary_TEs(batch_files, pan_genomes_dir, panTE_lib, output_dir, intact_ltr_paths, recover, log)

        gene_te_associations = output_dir + '/gene_te_associations.tsv'
        if len(gene_annotation_list) > 0:
            # 6. 根据 full_length_TE_gff 和 gene_gtf 获取插入到 gene 上、下游 1Kb、和内部的TE
            find_gene_relation_tes(batch_files, output_dir, recover, log)

        # 7. 根据RNA-seq数据识别由于LTR插入引起的差异表达基因
        if is_RNA_analyze:
            # 7.1 根据RNA-seq数据比对到基因组生成bam文件
            RNA_seq_dir = project_dir + '/RNA_seq'
            new_batch_files = generate_bam_for_RNA_seq(batch_files, threads, recover, RNA_seq_dir, log)

            # 7.2 调用 featureCount 进行定量
            gene_express_dir = output_dir + '/gene_quantities'
            gene_express_table = quantitative_gene(new_batch_files, gene_express_dir, threads, recover, log)

            # 7.3 调用R语言脚本 detect_DE_genes_from_TEs.R 找到由于LTR插入引起的差异表达基因
            script_dir = project_dir + '/RNA_seq'
            detect_DE_genes_from_TEs_cmd = 'cd ' + output_dir + ' && Rscript ' + script_dir + '/detect_DE_genes_from_TEs.R ' \
                                           + gene_express_table + ' ' + gene_te_associations
            os.system(detect_DE_genes_from_TEs_cmd)


if __name__ == '__main__':
    main_pipeline()
