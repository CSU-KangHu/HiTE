import argparse
import os
import sys
import time
from multiprocessing import cpu_count

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

    # 0. 读取 genomes
    genome_paths = []
    genome_names = []
    RNA_seq_dict = {}
    gene_annotation_list = []
    with open(genome_list_path, 'r') as f_r:
        for line in f_r:
            parts = line.replace('\n', '').split('\t')
            genome_name = parts[0]
            gene_gtf = None
            if len(parts) >= 2:
                gene_name = parts[1]
                if gene_name != '':
                    gene_gtf = genes_dir + '/' + gene_name
                    if not os.path.isabs(gene_gtf):
                        gene_gtf = os.path.abspath(gene_gtf)
                    if not os.path.exists(gene_gtf):
                        log.logger.error('\nCannot find gene annotation path: ' + gene_gtf)
                        exit(-1)
                    gene_annotation_list.append(gene_gtf)

            if len(parts) == 3:
                is_PE = False
                raw_RNA = RNA_dir + '/' + parts[2]
                if not os.path.exists(raw_RNA):
                    log.logger.error('\nCannot find RNA-seq path: ' + raw_RNA)
                    exit(-1)
                RNA_seq_dict['raw_RNA'] = raw_RNA
                RNA_seq_dict['is_PE'] = is_PE
            elif len(parts) == 4:
                is_PE = True
                raw_RNA1 = RNA_dir + '/' + parts[2]
                if not os.path.exists(raw_RNA1):
                    log.logger.error('\nCannot find RNA-seq path: ' + raw_RNA1)
                    exit(-1)
                raw_RNA2 = RNA_dir + '/' + parts[3]
                if not os.path.exists(raw_RNA2):
                    log.logger.error('\nCannot find RNA-seq path: ' + raw_RNA2)
                    exit(-1)
                RNA_seq_dict['raw_RNA1'] = raw_RNA1
                RNA_seq_dict['raw_RNA2'] = raw_RNA2
                RNA_seq_dict['is_PE'] = is_PE

            cur_genome_path = pan_genomes_dir + '/' + genome_name
            TE_gff = output_dir + '/' + genome_name + '.gff'

            if not os.path.isabs(cur_genome_path):
                cur_genome_path = os.path.abspath(cur_genome_path)
            if not os.path.exists(cur_genome_path):
                log.logger.error('\nCannot find genome path: ' + cur_genome_path)
                exit(-1)
            genome_paths.append((genome_name, cur_genome_path, TE_gff, gene_gtf, RNA_seq_dict))
            genome_names.append(genome_name)


    # 将gene注释文件转为标准的gtf格式，并去除非编码基因
    script_dir = project_dir + '/RNA_seq'
    genome_paths = convertGeneAnnotation2GTF(genome_paths, script_dir, output_dir, log)

    # 生成一个 total_genome.fa 用来进行 BM_EDTA 和 BM_HiTE 评估
    total_genome = output_dir + '/total_genome.fa'
    new_ref_contigs = {}
    for genome_name, reference, TE_gff, gene_gtf, RNA_seq_dict in genome_paths:
        raw_name = genome_name.split('.')[0]
        ref_names, ref_contigs = read_fasta(reference)
        for name in ref_names:
            new_name = raw_name + '-' + name
            new_ref_contigs[new_name] = ref_contigs[name]
    store_fasta(new_ref_contigs, total_genome)

    intact_ltr_paths = []
    pan_terminal_tmp_lib = output_dir + '/pan_terminal.tmp.fa'
    pan_internal_tmp_lib = output_dir + '/pan_internal.tmp.fa'
    panTE_lib = output_dir + '/panTE.fa'
    resut_file = panTE_lib
    if not recover or not os.path.exists(resut_file):
        os.system('rm -f ' + pan_terminal_tmp_lib)
        os.system('rm -f ' + pan_internal_tmp_lib)
        # 3. 利用 panTE 注释每个基因组
        for genome_name, reference, TE_gff, gene_gtf, RNA_seq_dict in genome_paths:
            raw_name = genome_name.split('.')[0]
            HiTE_output_dir = output_dir + '/HiTE_' + str(raw_name)
            ltr_intact_list = HiTE_output_dir + '/intact_LTR.list'
            intact_ltr_paths.append(ltr_intact_list)
            confident_ltr_terminal = HiTE_output_dir + '/confident_ltr.terminal.fa'
            confident_ltr_internal = HiTE_output_dir + '/confident_ltr.internal.fa'
            confident_helitron = HiTE_output_dir + '/confident_helitron.fa'
            confident_non_ltr = HiTE_output_dir + '/confident_non_ltr.fa'
            confident_other = HiTE_output_dir + '/confident_other.fa'
            confident_tir = HiTE_output_dir + '/confident_tir.fa'
            confident_TE = HiTE_output_dir + '/confident_TE.cons.fa'

            check_files = []
            is_rerun = True
            if te_type == 'all':
                if file_exist(confident_ltr_terminal) \
                        and file_exist(confident_ltr_internal) \
                        and file_exist(confident_helitron) \
                        and file_exist(confident_non_ltr) \
                        and file_exist(confident_other) \
                        and file_exist(confident_tir) \
                        and file_exist(confident_TE):
                    is_rerun = False
                    check_files.append(confident_ltr_terminal)
                    check_files.append(confident_ltr_internal)
                    check_files.append(confident_helitron)
                    check_files.append(confident_non_ltr)
                    check_files.append(confident_other)
                    check_files.append(confident_tir)
                    check_files.append(confident_TE)
            elif te_type == 'ltr':
                if file_exist(confident_ltr_terminal) \
                        and file_exist(confident_ltr_internal):
                    is_rerun = False
                    check_files.append(confident_ltr_terminal)
                    check_files.append(confident_ltr_internal)
            elif te_type == 'tir':
                if file_exist(confident_tir):
                    is_rerun = False
                    check_files.append(confident_tir)
            elif te_type == 'helitron':
                if file_exist(confident_tir):
                    is_rerun = False
                    check_files.append(confident_helitron)
            elif te_type == 'non-ltr':
                if file_exist(confident_other) \
                        and file_exist(confident_non_ltr):
                    is_rerun = False
                    check_files.append(confident_non_ltr)
                    check_files.append(confident_other)

            if not recover or is_rerun:
                if os.path.exists(reference):
                    # Step 2. 利用 HiTE 从泛基因组中检测TEs
                    HiTE_command = 'python ' + project_dir + '/main.py --genome ' + reference + ' --outdir ' + \
                                   HiTE_output_dir + ' --thread ' + str(threads) + ' --annotate 0 --te_type ' + \
                                   str(te_type) + ' --miu ' + str(miu) + ' --is_output_LTR_lib 0 ' + ' --debug ' + str(debug) + \
                                   ' --recover ' + str(recover)
                    log.logger.debug(HiTE_command)
                    starttime = time.time()
                    os.system(HiTE_command)
                    endtime = time.time()
                    dtime = endtime - starttime
                    log.logger.info("Running time of step2: %.4s m" % (dtime / 60))
                else:
                    log.logger.error('Cannot find Genome: ' + reference)
            else:
                for check_file in check_files:
                    log.logger.info(check_file + ' exists, skip...')


            # 对 library 加上前缀
            confident_ltr_terminal = lib_add_prefix(confident_ltr_terminal, raw_name)
            confident_ltr_internal = lib_add_prefix(confident_ltr_internal, raw_name)
            confident_TE = lib_add_prefix(confident_TE, raw_name)

            os.system('cat ' + confident_ltr_terminal + ' >> ' + pan_terminal_tmp_lib)
            os.system('cat ' + confident_ltr_internal + ' >> ' + pan_internal_tmp_lib)
            os.system('cat ' + confident_TE + ' >> ' + pan_terminal_tmp_lib)

            # # 删除临时目录
            # if not debug:
            #     os.system('rm -rf ' + HiTE_output_dir)

        # 对所有基因组的LTR终端和内部序列合集去冗余
        pan_terminal_tmp_lib_cons = pan_terminal_tmp_lib + '.cons'
        terminal_coverage_threshold = 0.95
        # Step7. Remove redundancy from the LTR terminal results.
        starttime = time.time()
        log.logger.info('Start step6: Remove LTR terminal redundancy')
        type = 'terminal'
        # 对于 terminal 而言，其内部不太可能出现太大的插入删除变异，因此我们只是利用已识别的LTR终端序列去冗余，不再获取拷贝支持
        deredundant_for_LTR_v5(pan_terminal_tmp_lib, output_dir, threads, type, terminal_coverage_threshold,
                               debug=0)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step6: %.8s s" % (dtime))


        pan_internal_tmp_lib_cons = pan_internal_tmp_lib + '.cons'
        internal_coverage_threshold = 0.8
        # Step7. Remove redundancy from the LTR internal results.
        starttime = time.time()
        log.logger.info('Start step6: Remove LTR internal redundancy')
        type = 'internal'
        # 对于 internal 而言，其内部可能出现较大的插入删除变异，因此我们要求比较宽松的0.8阈值
        deredundant_for_LTR_v5(pan_internal_tmp_lib, output_dir, threads, type, internal_coverage_threshold,
                               debug=0)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step6: %.8s s" % (dtime))

        # Step8. 生成一致性library
        os.system('cat ' + pan_terminal_tmp_lib_cons + ' ' + pan_internal_tmp_lib_cons + ' > ' + panTE_lib)
        # Reassign Inconsistent Classification Labels
        ReassignInconsistentLabels(panTE_lib)

    else:
        log.logger.info(resut_file + ' exists, skip...')
        for genome_name, reference, TE_gff, gene_gtf, RNA_seq_dict in genome_paths:
            raw_name = genome_name.split('.')[0]
            HiTE_output_dir = output_dir + '/HiTE_' + str(raw_name)
            ltr_intact_list = HiTE_output_dir + '/intact_LTR.list'
            intact_ltr_paths.append(ltr_intact_list)

    batch_files = []
    if not skip_analyze:
        # 3. 利用 panTE 注释每个基因组
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


        # Step 3. 根据注释好的 annotation 文件，进行常见TE分析
        log.logger.info('Start analysing using TE annotation files...')
        summary_TEs(batch_files, pan_genomes_dir, panTE_lib, output_dir, intact_ltr_paths, recover, log)

        gene_te_associations = output_dir + '/gene_te_associations.tsv'
        if len(gene_annotation_list) > 0:
            # Step 4. 根据 full_length_TE_gff 和 gene_gtf 获取插入到 gene 上、下游 1Kb、和内部的TE
            find_gene_relation_tes(batch_files, output_dir, recover, log)

        if len(RNA_seq_dict) > 0:
            # Step 5. 根据RNA-seq数据比对到基因组生成bam文件
            RNA_seq_dir = project_dir + '/RNA_seq'
            new_batch_files = generate_bam_for_RNA_seq(batch_files, threads, recover, RNA_seq_dir, log)

            # Step 6. 调用 featureCount 进行定量
            gene_express_dir = output_dir + '/gene_quantities'
            gene_express_table = quantitative_gene(new_batch_files, gene_express_dir, threads, recover, log)

            # step 7. 调用R语言脚本 detect_DE_genes_from_TEs.R 找到由于LTR插入引起的差异表达基因
            script_dir = project_dir + '/RNA_seq'
            detect_DE_genes_from_TEs_cmd = 'cd ' + output_dir + ' && Rscript ' + script_dir + '/detect_DE_genes_from_TEs.R ' \
                                           + gene_express_table + ' ' + gene_te_associations
            os.system(detect_DE_genes_from_TEs_cmd)


if __name__ == '__main__':
    main_pipeline()
