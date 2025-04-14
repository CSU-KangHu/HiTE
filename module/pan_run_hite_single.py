#!/usr/bin/env python
import argparse
import os
import shutil
import time
import uuid

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import file_exist, lib_add_prefix, Logger, store_fasta, copy_files, create_or_clear_directory, \
    clean_old_tmp_files_by_dir


def for_test(genome_name, reference, output_dir, threads, te_type, miu, debug, log):
    """对单个基因组运行 HiTE"""
    raw_name = genome_name.split('.')[0]
    HiTE_output_dir = output_dir

    # 定义所需的文件路径
    ltr_intact_list = os.path.join(HiTE_output_dir, "intact_LTR.list")
    confident_ltr_terminal = os.path.join(HiTE_output_dir, "confident_ltr.terminal.fa")
    confident_ltr_internal = os.path.join(HiTE_output_dir, "confident_ltr.internal.fa")
    confident_ltr_cut = os.path.join(HiTE_output_dir, "confident_ltr_cut.fa")
    intact_ltr = os.path.join(HiTE_output_dir, "intact_LTR.fa")
    intact_ltr_classified = os.path.join(HiTE_output_dir, "intact_LTR.fa.classified")
    confident_helitron = os.path.join(HiTE_output_dir, "confident_helitron.fa")
    confident_non_ltr = os.path.join(HiTE_output_dir, "confident_non_ltr.fa")
    confident_other = os.path.join(HiTE_output_dir, "confident_other.fa")
    confident_tir = os.path.join(HiTE_output_dir, "confident_tir.fa")
    confident_TE = os.path.join(HiTE_output_dir, "confident_TE.cons.fa")

    file_paths = [ltr_intact_list, confident_ltr_terminal, confident_ltr_internal, confident_ltr_cut, intact_ltr,
                  intact_ltr_classified, confident_helitron, confident_non_ltr, confident_other, confident_tir,
                  confident_TE]
    for file_path in file_paths:
        os.system('echo aaaa > ' + file_path)



def run_hite_for_genome(genome_name, reference, output_dir, threads, te_type, miu, shared_prev_TE, debug, log):
    """对单个基因组运行 HiTE"""
    raw_name = genome_name.split('.')[0]
    HiTE_output_dir = output_dir

    # 定义所需的文件路径
    ltr_intact_list = os.path.join(HiTE_output_dir, "intact_LTR.list")
    confident_ltr_terminal = os.path.join(HiTE_output_dir, "confident_ltr.terminal.fa")
    confident_ltr_internal = os.path.join(HiTE_output_dir, "confident_ltr.internal.fa")
    confident_ltr_cut_path = os.path.join(HiTE_output_dir, "confident_ltr_cut.fa")

    confident_helitron = os.path.join(HiTE_output_dir, "confident_helitron.fa")
    confident_non_ltr = os.path.join(HiTE_output_dir, "confident_non_ltr.fa")
    confident_other = os.path.join(HiTE_output_dir, "confident_other.fa")
    confident_tir = os.path.join(HiTE_output_dir, "confident_tir.fa")
    confident_TE = os.path.join(HiTE_output_dir, "confident_TE.cons.fa")

    add_prefix_files = [confident_ltr_terminal, confident_ltr_internal, confident_ltr_cut_path,
                        confident_helitron, confident_non_ltr, confident_other, confident_tir, confident_TE]

    # 运行 HiTE
    if file_exist(reference):
        HiTE_command = (
            f"python {project_dir}/main.py --genome {reference} --out_dir {HiTE_output_dir} "
            f"--thread {threads} --annotate 0 --te_type {te_type} --miu {miu} --is_output_LTR_lib 1 --recover 1  "
            f"--debug {debug} --work_dir {HiTE_output_dir}"
        )
        if shared_prev_TE is not None:
            HiTE_command += f" --shared_prev_TE {shared_prev_TE}"
        log.logger.info(f"Executing: {HiTE_command}")
        start_time = time.time()
        os.system(HiTE_command)
        end_time = time.time()
        log.logger.info(f"Running time for {genome_name}: {(end_time - start_time) / 60:.2f} minutes")
        # 为文件加前缀
        for cur_file in add_prefix_files:
            if file_exist(cur_file):
                lib_add_prefix(cur_file, raw_name)

    else:
        log.logger.error(f"Cannot find genome: {reference}")




def main(genome_name, reference, output_dir, threads, te_type, miu, shared_prev_TE, debug, log):
    """主函数"""
    # 运行单个基因组的 HiTE 检测
    run_hite_for_genome(genome_name, reference, output_dir, threads, te_type, miu, shared_prev_TE, debug, log)

    # for_test(genome_name, reference, output_dir, threads, te_type, miu, debug, log)

if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE run single genome.")
    parser.add_argument("--genome_name", type=str, help="Name of the genome.")
    parser.add_argument("--reference", type=str, help="Path to the reference file.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--te_type", type=str, help="Type of transposable element (TE).")
    parser.add_argument("--miu", type=float, help="Parameter miu for the process.")
    parser.add_argument("--debug", type=int, help="Enable or disable debug mode (True/False).")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(), help="Output directory (default: current working directory).")
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")
    parser.add_argument('--shared_prev_TE', metavar='shared_prev_TE', help='The path of shared previous TEs')

    # 解析参数
    args = parser.parse_args()
    genome_name = args.genome_name
    reference = args.reference
    threads = args.threads
    te_type = args.te_type
    miu = args.miu
    debug = args.debug
    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)
    shared_prev_TE = args.shared_prev_TE

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'pan_run_hite_single_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)

        main(genome_name, reference, temp_dir, threads, te_type, miu, shared_prev_TE, debug, log)

        # 计算完之后将结果拷贝回输出目录
        copy_files(temp_dir, output_dir)

    except Exception as e:
        # 如果出现异常，打印错误信息并删除临时目录
        print(f"An error occurred: {e}")
        # if os.path.exists(temp_dir):
        #     shutil.rmtree(temp_dir)
        raise  # 重新抛出异常，以便上层代码可以处理

    else:
        # 如果没有异常，删除临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
