import argparse
import json
import os
import sys
import time
from multiprocessing import cpu_count

current_folder = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.join(current_folder, ".")

from module.Util import Logger, file_exist


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
    parser.add_argument('--genes_dir', required=False, metavar='genes_dir', help='A directory containing the gene annotation files, gff or gtf format.')
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
    genome_metadata = os.path.join(output_dir, 'genome_metadata.json')
    result_file = genome_metadata
    if not recover or not file_exist(result_file):
        starttime = time.time()
        pan_preprocess_cmd = 'pan_preprocess_genomes.py ' + genome_list_path + ' ' + genes_dir \
                             + ' ' + RNA_dir + ' ' + pan_genomes_dir + ' ' + output_dir
        log.logger.info(pan_preprocess_cmd)
        os.system(pan_preprocess_cmd)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step1: %.8s s" % (dtime))
    else:
        log.logger.info(result_file + ' exists, skip...')

    # Load the metadata
    with open(genome_metadata, 'r') as f:
        genome_data = json.load(f)
    genome_info_list = genome_data["genome_info"]

    # 2. 遍历每个基因组，并调用HiTE进行TE识别
    intact_ltr_paths = []
    pan_terminal_tmp_lib = os.path.join(output_dir, 'pan_terminal.tmp.fa')
    pan_internal_tmp_lib = os.path.join(output_dir, 'pan_internal.tmp.fa')
    panTE_lib = os.path.join(output_dir, 'panTE.fa')

    result_file = panTE_lib
    if not recover or not file_exist(result_file):
        starttime = time.time()
        for genome_info in genome_info_list:
            genome_name = genome_info["genome_name"]
            raw_name = genome_info["raw_name"]
            reference = genome_info["reference"]
            HiTE_output_dir = os.path.join(output_dir, f'HiTE_{raw_name}')
            single_result_path = os.path.join(HiTE_output_dir, f"{genome_name}_hite_result.json")
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
        log.logger.info(result_file + ' exists, skip...')
        for genome_info in genome_info_list:
            raw_name = genome_info["raw_name"]
            HiTE_output_dir = os.path.join(output_dir, f'HiTE_{raw_name}')
            ltr_intact_list = HiTE_output_dir + '/intact_LTR.list'
            intact_ltr_paths.append(ltr_intact_list)

    intact_ltr_paths_file = os.path.join(output_dir, 'intact_ltr_paths.txt')
    with open(intact_ltr_paths_file, 'w') as f:
        json.dump(intact_ltr_paths, f, indent=4)

    # 4. 利用 panTE 注释每个基因组
    for genome_info in genome_info_list:
        genome_name = genome_info["genome_name"]
        reference = genome_info["reference"]
        TE_gff = genome_info["TE_gff"]
        pan_annotate_genome_cmd = 'pan_annotate_genome.py ' + TE_gff + ' ' + output_dir + ' ' + str(threads) \
                                  + ' ' + panTE_lib + ' ' + reference + ' ' + genome_name + ' ' + str(recover)
        log.logger.info(pan_annotate_genome_cmd)
        os.system(pan_annotate_genome_cmd)

    if not skip_analyze:
        # 5. 根据注释好的 annotation 文件，进行常见TE分析
        pan_summary_TEs_cmd = 'pan_summary_TEs.py ' + genome_metadata + ' ' + pan_genomes_dir + ' ' + panTE_lib \
                                  + ' ' + output_dir + ' ' + intact_ltr_paths_file + ' ' + str(recover)
        log.logger.info(pan_summary_TEs_cmd)
        os.system(pan_summary_TEs_cmd)

        # 6. 根据 full_length_TE_gff 和 gene_gtf 获取插入到 gene 上、下游 1Kb、和内部的TE
        gene_te_associations = output_dir + '/gene_te_associations.tsv'
        pan_gene_te_relation_cmd = 'pan_gene_te_relation.py ' + genome_metadata + ' ' + output_dir + ' ' + str(recover)
        log.logger.info(pan_gene_te_relation_cmd)
        os.system(pan_gene_te_relation_cmd)


        # 7. 根据RNA-seq数据识别由于LTR插入引起的差异表达基因
        pan_detect_de_genes_cmd = 'pan_detect_de_genes.py ' + genome_metadata + ' ' + str(threads) + ' ' + str(recover) \
                                  + ' ' + output_dir  + ' ' + gene_te_associations
        log.logger.info(pan_detect_de_genes_cmd)
        os.system(pan_detect_de_genes_cmd)


if __name__ == '__main__':
    main_pipeline()
