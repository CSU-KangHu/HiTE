import os
import sys
import time
import json
import argparse
from Util import file_exist, lib_add_prefix, Logger

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

def run_hite_for_genome(genome_name, reference, output_dir, threads, te_type, miu, debug, recover, log):
    """对单个基因组运行 HiTE"""
    raw_name = genome_name.split('.')[0]
    HiTE_output_dir = os.path.join(output_dir, f"HiTE_{raw_name}")

    # 定义所需的文件路径
    ltr_intact_list = os.path.join(HiTE_output_dir, "intact_LTR.list")
    confident_ltr_terminal = os.path.join(HiTE_output_dir, "confident_ltr.terminal.fa")
    confident_ltr_internal = os.path.join(HiTE_output_dir, "confident_ltr.internal.fa")
    confident_helitron = os.path.join(HiTE_output_dir, "confident_helitron.fa")
    confident_non_ltr = os.path.join(HiTE_output_dir, "confident_non_ltr.fa")
    confident_other = os.path.join(HiTE_output_dir, "confident_other.fa")
    confident_tir = os.path.join(HiTE_output_dir, "confident_tir.fa")
    confident_TE = os.path.join(HiTE_output_dir, "confident_TE.cons.fa")

    # 初始化文件检查列表
    check_files = []
    if te_type == "all":
        check_files = [
            confident_ltr_terminal, confident_ltr_internal, confident_helitron,
            confident_non_ltr, confident_other, confident_tir, confident_TE
        ]
    elif te_type == "ltr":
        check_files = [confident_ltr_terminal, confident_ltr_internal]
    elif te_type == "tir":
        check_files = [confident_tir]
    elif te_type == "helitron":
        check_files = [confident_helitron]
    elif te_type == "non-ltr":
        check_files = [confident_non_ltr, confident_other]

    # 判断是否需要重新运行 HiTE
    is_rerun = not all(file_exist(f) for f in check_files)

    # 运行 HiTE
    if not recover or is_rerun:
        if file_exist(reference):
            HiTE_command = (
                f"python {project_dir}/main.py --genome {reference} --outdir {HiTE_output_dir} "
                f"--thread {threads} --annotate 0 --te_type {te_type} --miu {miu} --is_output_LTR_lib 0 "
                f"--debug {debug} --recover {recover}"
            )
            log.logger.info(f"Executing: {HiTE_command}")
            start_time = time.time()
            os.system(HiTE_command)
            end_time = time.time()
            log.logger.info(f"Running time for {genome_name}: {(end_time - start_time) / 60:.2f} minutes")

            # 为文件加前缀
            if file_exist(confident_ltr_terminal):
                confident_ltr_terminal = lib_add_prefix(confident_ltr_terminal, raw_name)
            if file_exist(confident_ltr_internal):
                confident_ltr_internal = lib_add_prefix(confident_ltr_internal, raw_name)
            if file_exist(confident_TE):
                confident_TE = lib_add_prefix(confident_TE, raw_name)
        else:
            log.logger.error(f"Cannot find genome: {reference}")
    else:
        for check_file in check_files:
            log.logger.info(f"{check_file} exists, skipping HiTE run...")

    # 返回所有生成文件路径
    return {
        "genome_name": genome_name,
        "ltr_intact_list": ltr_intact_list if file_exist(ltr_intact_list) else None,
        "confident_ltr_terminal": confident_ltr_terminal if file_exist(confident_ltr_terminal) else None,
        "confident_ltr_internal": confident_ltr_internal if file_exist(confident_ltr_internal) else None,
        "confident_helitron": confident_helitron if file_exist(confident_helitron) else None,
        "confident_non_ltr": confident_non_ltr if file_exist(confident_non_ltr) else None,
        "confident_other": confident_other if file_exist(confident_other) else None,
        "confident_tir": confident_tir if file_exist(confident_tir) else None,
        "confident_TE": confident_TE if file_exist(confident_TE) else None,
    }


def main(genome_name, reference, output_dir, threads, te_type, miu, debug, recover, log):
    """主函数"""
    # 运行单个基因组的 HiTE 检测
    result = run_hite_for_genome(genome_name, reference, output_dir, threads, te_type, miu, debug, recover, log)

    # 保存结果到 JSON
    result_path = os.path.join(output_dir, f"{genome_name}_hite_result.json")
    with open(result_path, 'w') as f:
        json.dump(result, f, indent=4)
    log.logger.info(f"Result saved to: {result_path}")

if __name__ == "__main__":
    # 接收命令行参数
    genome_name = sys.argv[1]
    reference = sys.argv[2]
    output = sys.argv[3]
    threads = sys.argv[4]
    te_type = sys.argv[5]
    miu = sys.argv[6]
    debug = sys.argv[7]
    recover = sys.argv[8]

    if output is None:
        output = os.getcwd()
    output = os.path.abspath(output)
    log = Logger(output + '/panHiTE.log', level='debug')
    main(genome_name, reference, output, threads, te_type, miu, debug, recover, log)
