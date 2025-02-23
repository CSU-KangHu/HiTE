import argparse
import os
import sys


current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from Util import Logger, rename_reference, file_exist

from configs import config

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run benchmarking...')
    parser.add_argument('--BM_RM2', metavar='BM_RM2',
                        help='Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false.')
    parser.add_argument('--BM_EDTA', metavar='BM_EDTA',
                        help='Whether to conduct benchmarking of EDTA, 1: true, 0: false.')
    parser.add_argument('--BM_HiTE', metavar='BM_HiTE',
                        help='Whether to conduct benchmarking of HiTE, 1: true, 0: false.')
    parser.add_argument('--coverage_threshold', metavar='coverage_threshold',
                        help='coverage threshold')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--TE_lib', metavar='TE_lib',
                        help='The path of TE library')
    parser.add_argument('-r', metavar='reference',
                        help='e.g., Input reference Path')
    parser.add_argument('--EDTA_home', metavar='EDTA_home',
                        help='When conducting benchmarking of EDTA, you will be asked to input EDTA home path.')
    parser.add_argument('--species', metavar='species',
                        help='Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--recover', metavar='is_recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')


    args = parser.parse_args()
    BM_RM2 = args.BM_RM2
    BM_EDTA = args.BM_EDTA
    BM_HiTE = args.BM_HiTE
    threads = int(args.t)
    TE_lib = args.TE_lib
    reference = args.r
    EDTA_home = args.EDTA_home
    species = args.species
    tmp_output_dir = args.tmp_output_dir
    coverage_threshold = args.coverage_threshold
    recover = int(args.recover)

    default_coverage_threshold = 0.95

    if coverage_threshold is not None:
        coverage_threshold = float(coverage_threshold)
    else:
        coverage_threshold = default_coverage_threshold

    project_dir = config.project_dir
    src_dir = project_dir + '/src'

    rm2_script = project_dir + '/tools/get_family_summary_paper.sh'
    rm2_script_improved = project_dir + '/tools/get_family_summary_paper_improved.sh'
    rm2_strict_script = project_dir + '/tools/get_family_summary_paper_0.99.sh'
    lib_module = project_dir + '/library'
    if species == "dmel":
        lib_path = lib_module + "/dmel.ltr.ref"
    elif species == "rice":
        lib_path = lib_module + "/rice.ltr.ref"
    elif species == "cb":
        lib_path = lib_module + "/cbrrep.ref"
    elif species == "zebrafish":
        lib_path = lib_module + "/zebrafish.ltr.ref"
    elif species == "maize":
        lib_path = lib_module + "/maize.ltr.ref"
    elif species == "ath":
        lib_path = lib_module + "/ath.ltr.ref"
    elif species == "chicken":
        lib_path = lib_module + "/chicken.ref"
    elif species == "zebrafinch":
        lib_path = lib_module + "/zebrafinch.ref"
    elif species == "mouse":
        lib_path = lib_module + "/mouse.ltr.ref"
    elif species == "human":
        lib_path = lib_module + "/human.ltr.ref"
    elif species == "Chlamydomonas_reinhardtii":
        lib_path = lib_module + "/Chlamydomonas_reinhardtii.ltr.ref"
    elif species == "Picea_abies":
        lib_path = lib_module + "/Picea_abies.ltr.ref"
    elif species == "Thalassiosira_pseudonana":
        lib_path = lib_module + "/Thalassiosira_pseudonana.ltr.ref"
    else:
        lib_path = lib_module + "/test.ref"

    curated_lib = os.path.abspath(lib_path)
    TE_lib = os.path.abspath(TE_lib)
    reference = os.path.abspath(reference)
    rm2_script = os.path.abspath(rm2_script)
    rm2_strict_script = os.path.abspath(rm2_strict_script)

    if coverage_threshold == 0.99:
        rm2_script = rm2_strict_script

    log = Logger(tmp_output_dir+'/benchmarking.log', level='debug')

    rm2_test_dir = tmp_output_dir + '/rm2_test'
    rm2_out = tmp_output_dir + '/BM_RM2.log'
    edta_out = tmp_output_dir + '/BM_EDTA.log'
    hite_out = tmp_output_dir + '/BM_HiTE.log'

    os.system('rm -rf ' + rm2_test_dir)

    is_BM_RM2 = False
    is_BM_EDTA = False
    is_BM_HiTE = False
    if BM_RM2 is not None and int(BM_RM2) == 1:
        is_BM_RM2 = True
    if (BM_EDTA is not None and int(BM_EDTA) == 1):
        is_BM_EDTA = True
    if (BM_HiTE is not None and int(BM_HiTE) == 1):
        is_BM_HiTE = True

    # rename reference
    ref_rename_path = tmp_output_dir + '/genome.rename.fa'
    chr_name_map = tmp_output_dir + '/chr_name.map'
    rename_reference(reference, ref_rename_path, chr_name_map)
    reference = ref_rename_path

    # 2.annotate the genome
    if is_BM_RM2:
        rm_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -lib ' + curated_lib + ' -nolow -pa ' + str(threads) + ' ' + TE_lib
        log.logger.debug(rm_command)
        os.system(rm_command)
        if not os.path.exists(rm2_test_dir):
            os.makedirs(rm2_test_dir)
        result_command = 'cd ' + rm2_test_dir + ' && sh ' + rm2_script + ' ' + TE_lib + '.out ' + ' >> ' + rm2_out
        log.logger.debug(result_command)
        os.system(result_command)

        # if os.path.exists(rm2_test_dir):
        #     os.system('rm -rf ' + rm2_test_dir)
        # if not os.path.exists(rm2_test_dir):
        #     os.makedirs(rm2_test_dir)
        # result_command = 'cd ' + rm2_test_dir + ' && sh ' + rm2_script_improved + ' ' + TE_lib + '.out ' + curated_lib + ' ' + TE_lib + ' >> ' + rm2_out
        # log.logger.debug(result_command)
        # os.system(result_command)
    else:
        log.logger.debug('Skip benchmarking of RepeatModeler2')

    if is_BM_HiTE:
        divergence_threshold = (1 - coverage_threshold) * 100

        bm_hite_command = 'cd ' + src_dir + ' && python ' + src_dir + '/lib_evaluation.py -g ' \
                          + reference + ' --standard_lib ' + lib_path \
                          + ' --test_lib ' + TE_lib + ' --work_dir ' + tmp_output_dir \
                          + ' --coverage_threshold ' + str(coverage_threshold) + ' --thread ' + str(threads) \
                          + ' ' + ' --cat Total ' + ' --is_full_length 1 ' + ' >> ' + hite_out
        log.logger.debug(bm_hite_command)
        os.system(bm_hite_command)
    else:
        log.logger.debug('Skip benchmarking of HiTE')

    if is_BM_EDTA:
        repbase_out = tmp_output_dir + '/repbase.edta.out'
        result_file = repbase_out
        if not recover or not file_exist(result_file):
            RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                                   + ' -q -no_is -norna -nolow -div 40 -gff -lib ' + curated_lib + ' -cutoff 225 ' \
                                   + reference
            log.logger.debug(RepeatMasker_command)
            os.system(RepeatMasker_command)

            mv_file_command = 'mv ' + reference + '.out ' + repbase_out
            log.logger.debug(mv_file_command)
            os.system(mv_file_command)
        else:
            log.logger.info(result_file + ' exists, skip...')

        test_out = tmp_output_dir + '/HiTE.edta.out'
        result_file = test_out
        if not recover or not file_exist(result_file):
            RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                                   + ' -q -no_is -norna -nolow -div 40 -gff -lib ' + TE_lib + ' -cutoff 225 ' \
                                   + reference
            log.logger.debug(RepeatMasker_command)
            os.system(RepeatMasker_command)

            mv_file_command = 'mv ' + reference + '.out ' + test_out
            log.logger.debug(mv_file_command)
            os.system(mv_file_command)
        else:
            log.logger.info(result_file + ' exists, skip...')

        # remove old report
        os.system('rm -f ' + tmp_output_dir + '/HiTE.edta.out.*.report')
        bm_edta_command = 'cd ' + tmp_output_dir + ' && perl ' + EDTA_home + '/lib-test.pl -genome ' \
                          + reference + ' -std ' + repbase_out + ' -tst ' + test_out + ' -cat Total'
        log.logger.debug(bm_edta_command)
        os.system(bm_edta_command)

        mv_file_command = 'mv ' + tmp_output_dir + '/HiTE.edta.out.*.report ' + edta_out
        log.logger.debug(mv_file_command)
        os.system(mv_file_command)
    else:
        log.logger.debug('Skip benchmarking of EDTA')




