import argparse
import os
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import Logger



if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE benchmarking...')
    parser.add_argument('--BM_RM2', metavar='BM_RM2',
                        help='Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false.')
    parser.add_argument('--BM_EDTA', metavar='BM_EDTA',
                        help='Whether to conduct benchmarking of EDTA, 1: true, 0: false.')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--lib_module', metavar='lib_module',
                        help='Please enter the directory for library. Use an absolute path.')
    parser.add_argument('--TE_lib', metavar='TE_lib',
                        help='The path of TE library')
    parser.add_argument('--rm2_script', metavar='rm2_script',
                        help='The path of BM_RM2 script')
    parser.add_argument('-r', metavar='reference',
                        help='e.g., Input reference Path')
    parser.add_argument('--EDTA_home', metavar='EDTA_home',
                        help='When conducting benchmarking of EDTA, you will be asked to input EDTA home path.')
    parser.add_argument('--species', metavar='species',
                        help='Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')


    args = parser.parse_args()
    BM_RM2 = args.BM_RM2
    BM_EDTA = args.BM_EDTA
    threads = int(args.t)
    lib_module = args.lib_module
    TE_lib = args.TE_lib
    rm2_script = args.rm2_script
    reference = args.r
    EDTA_home = args.EDTA_home
    species = args.species
    tmp_output_dir = args.tmp_output_dir

    if species == "dmel":
        lib_path = lib_module + "/drorep.ref"
    elif species == "rice":
        lib_path = lib_module + "/oryrep.ref"
    elif species == "cb":
        lib_path = lib_module + "/cbrrep.ref"
    elif species == "zebrafish":
        lib_path = lib_module + "/zebrep.ref"
    elif species == "maize":
        lib_path = lib_module + "/maize.ref"
    elif species == "ath":
        lib_path = lib_module + "/athrep.ref"
    else:
        lib_path = lib_module + "/test.ref"

    curated_lib = os.path.abspath(lib_path)
    TE_lib = os.path.abspath(TE_lib)
    reference = os.path.abspath(reference)
    rm2_script = os.path.abspath(rm2_script)

    log = Logger(tmp_output_dir+'/benchmarking.log', level='debug')

    rm2_test_dir = tmp_output_dir + '/rm2_test'
    res_out = tmp_output_dir + '/BM_RM2.log'

    os.system('rm -rf ' + rm2_test_dir)


    # 2.annotate the genome
    if BM_RM2 is not None and int(BM_RM2) == 1:
        rm_command = 'RepeatMasker -lib ' + curated_lib + ' -nolow -pa ' + str(threads) + ' ' + TE_lib
        log.logger.debug(rm_command)
        os.system(rm_command)
        if not os.path.exists(rm2_test_dir):
            os.makedirs(rm2_test_dir)
        result_command = 'cd ' + rm2_test_dir + ' && sh ' + rm2_script + ' ' + TE_lib + '.out > ' + res_out
        log.logger.debug(result_command)
        os.system(result_command)
    else:
        log.logger.debug('Skip benchmarking of RepeatModeler2')

    if BM_EDTA is not None and int(BM_EDTA) == 1:
        RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                               + ' -q -no_is -norna -nolow -div 40 -gff -lib ' + curated_lib + ' -cutoff 225 ' \
                               + reference
        log.logger.debug(RepeatMasker_command)
        os.system(RepeatMasker_command)

        mv_file_command = 'mv ' + reference + '.out ' + tmp_output_dir + '/repbase.out'
        log.logger.debug(mv_file_command)
        os.system(mv_file_command)

        RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                               + ' -q -no_is -norna -nolow -div 40 -gff -lib ' + TE_lib + ' -cutoff 225 ' \
                               + reference
        log.logger.debug(RepeatMasker_command)
        os.system(RepeatMasker_command)

        mv_file_command = 'mv ' + reference + '.out ' + tmp_output_dir + '/HiTE.out'
        log.logger.debug(mv_file_command)
        os.system(mv_file_command)

        bm_edta_command = 'cd ' + tmp_output_dir + ' && perl ' + EDTA_home + '/lib-test.pl -genome ' + reference + ' -std ' + tmp_output_dir \
                          + '/repbase.out' + ' -tst ' + tmp_output_dir + '/HiTE.out' + ' -cat Total'
        log.logger.debug(bm_edta_command)
        os.system(bm_edta_command)

        mv_file_command = 'mv ' + tmp_output_dir + '/HiTE.out.*.report ' + tmp_output_dir + '/BM_EDTA.log'
        log.logger.debug(mv_file_command)
        os.system(mv_file_command)
    else:
        log.logger.debug('Skip benchmarking of EDTA')


