#-- coding: UTF-8 --
import argparse
import codecs
import datetime
import json
import multiprocessing
import os
import shutil
import time

def rename_header(fasta_path):
    contigs = {}
    contigNames = []
    node_index = 0
    with open(fasta_path, "r") as f_r:
        contigName = ''
        contigseq = ''
        for line in f_r:
            if line.startswith('>'):
                if contigName != '' and contigseq != '':
                    contigs[contigName] = contigseq
                    contigNames.append(contigName)
                contigName = "Node_" + str(node_index)
                node_index += 1
                contigseq = ''
            else:
                contigseq += line.strip().upper()
        contigs[contigName] = contigseq
        contigNames.append(contigName)
    f_r.close()
    return contigNames, contigs

def read_fasta(fasta_path):
    contigs = {}
    contigNames = []
    with open(fasta_path, "r") as f_r:
        contigName = ''
        contigseq = ''
        for line in f_r:
            if line.startswith('>'):
                if contigName != '' and contigseq != '':
                    contigs[contigName] = contigseq
                    contigNames.append(contigName)
                contigName = line.strip()[1:].split(' ')[0]
                contigseq = ''
            else:
                contigseq += line.strip().upper()
        contigs[contigName] = contigseq
        contigNames.append(contigName)
    f_r.close()
    return contigNames, contigs

def read_classified_fasta(fasta_path):
    unclass_class_names_map = {}
    contigs = {}
    contigNames = []
    with open(fasta_path, "r") as f_r:
        contigName = ''
        contigseq = ''
        node_name = ''
        class_name = ''
        node_info = ''
        for line in f_r:
            if line.startswith('>'):
                if contigName != '' and contigseq != '':
                    contigs[contigName] = contigseq
                    contigNames.append(contigName)
                contigName = line.strip()[1:]

                classified_index = contigName.find('#')
                blank_index = contigName.find(' ')
                node_name = contigName[0:classified_index]
                if blank_index != -1:
                    class_name = contigName[classified_index+1: blank_index]
                    node_info = ' ' + contigName[blank_index+1:]
                else:
                    class_name = contigName[classified_index + 1: ]
                    node_info = ''
                unclass_name = node_name + node_info
                unclass_class_names_map[unclass_name] = contigName

                contigseq = ''
            else:
                contigseq += line.strip().upper()
        contigs[contigName] = contigseq
        contigNames.append(contigName)
    f_r.close()
    return (unclass_class_names_map, contigNames, contigs)

def modify_REPCLASS_conf(repclass_conf_path, repclass_conf):
    new_str = ''
    with open(repclass_conf_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('$'):
                key = line.split('=')[0].strip()[1:]
                if repclass_conf.__contains__(key):
                    value = repclass_conf[key]
                    line = '$'+key+'		=	\"'+value+'\";'
            new_str += line + '\n'
    f_r.close()

    with open(repclass_conf_path, 'w') as f_save:
        f_save.write(new_str)
        f_save.flush()
    f_save.close()

def set_REPCLASS_conf(cur_user_config_dir, cur_sample_name, genome, cur_consensus, cur_output, partition_index):
    create_dirs = []
    create_dirs.append(cur_user_config_dir)
    create_dirs.append(cur_output)
    for create_dir in create_dirs:
        if not os.path.exists(create_dir):
            os.makedirs(create_dir)

    # 1. copy original config file
    original_user_config_dir = os.getcwd() + "/third-party/REPCLASS-master/1.0.1/Myconf"
    os.system('cp -r ' + original_user_config_dir + ' ' + cur_user_config_dir)


    # 2. set userconfigure.conf
    userconfigure_conf = {}
    userconfigure_conf_path = cur_user_config_dir + '/Myconf/userconfigure.conf'
    repclass_conf_path = os.getcwd() + "/third-party/REPCLASS-master/1.0.1/conf/repclass.conf"
    JOB_NAME = cur_sample_name
    DATA = cur_output
    TEMPCHOICE = cur_user_config_dir + "/Myconf/mytempchoice"
    (GENOME_LOCATION, GENOME_FILE) = os.path.split(genome)
    TE_SEQUENCE = cur_consensus
    REPCLASS_CONF = repclass_conf_path

    userconfigure_conf['JOB_NAME'] = JOB_NAME
    userconfigure_conf['DATA'] = DATA
    userconfigure_conf['TEMPCHOICE'] = TEMPCHOICE
    userconfigure_conf['GENOME_LOCATION'] = GENOME_LOCATION
    userconfigure_conf['GENOME_FILE'] = GENOME_FILE
    userconfigure_conf['TE_SEQUENCE'] = TE_SEQUENCE
    userconfigure_conf['REPCLASS_CONF'] = REPCLASS_CONF
    modify_REPCLASS_conf(userconfigure_conf_path, userconfigure_conf)


#seq_item -> (header, sequence)
def PET(seq_item, partitions):
    # sort contigs by length
    original = seq_item
    original = sorted(original, key=lambda x: len(x[1]), reverse=True)
    return divided_array(original, partitions)


def divided_array(original_array, partitions):
    final_partitions = [[] for _ in range(partitions)]
    node_index = 0

    read_from_start = True
    read_from_end = False
    i = 0
    j = len(original_array) - 1
    while i <= j:
        # read from file start
        if read_from_start:
            final_partitions[node_index % partitions].append(original_array[i])
            i += 1
        if read_from_end:
            final_partitions[node_index % partitions].append(original_array[j])
            j -= 1
        node_index += 1
        if node_index % partitions == 0:
            # reverse
            read_from_end = bool(1 - read_from_end)
            read_from_start = bool(1 - read_from_start)
    return final_partitions


def store2file(data_partition, cur_consensus_path):
    if len(data_partition) > 0:
        with open(cur_consensus_path, 'w') as f_save:
            for item in data_partition:
                f_save.write('>'+item[0]+'\n'+item[1]+'\n')
        f_save.close()

def getREPCLASS_classified(REPCLASS_result_path):
    marked_header = []
    if os.path.exists(REPCLASS_result_path):
        with open(REPCLASS_result_path, 'r') as f_r:
            for index, line in enumerate(f_r):
                line = line.replace('\n', '').strip()
                if index == 0 or line == '':
                    continue
                parts = line.split('\t')
                class_name = 'Unknown'
                node_name = ''
                if len(parts) == 1:
                    node_name = parts[0].strip()
                elif len(parts) > 1:
                    node_name = parts[0].strip()
                    Class1 = parts[7].strip()
                    if Class1 == 'S':
                        super_family = parts[9].strip()
                        family = parts[15].strip()
                        class_name = super_family
                        if family != '':
                            class_name += '/' + family
                    elif Class1 == 'T':
                        Class2 = parts[10].strip()
                        if Class2 == 'S':
                            super_family = parts[12].strip()
                            family = parts[18].strip()
                            class_name = super_family
                            if family != '':
                                class_name += '/' + family
                        else:
                            super_family = parts[9].strip()
                            family = parts[15].strip()
                            class_name = super_family
                            if family != '':
                                class_name += '/' + family
                # mark each sequence
                marked_header.append(node_name + '#' + class_name)
        f_r.close()
    return marked_header


def merge_classified(marked_header, homology_classified_path, homology_unknown_path, final_classified_path):
    classified_contigNames, classified_contigs = read_fasta(homology_classified_path)
    unknown_contigNames, unknown_contigs = read_fasta(homology_unknown_path)
    with open(final_classified_path, 'w') as f_save:
        for name in classified_contigNames:
            f_save.write('>' + name + '\n' + classified_contigs[name] + '\n')
        for name in marked_header:
            header = name.split('#')[0]
            if unknown_contigs.__contains__(header):
                f_save.write('>' + name + '\n' + unknown_contigs[header] + '\n')
    f_save.close()


def run_classification(cur_consensus_path, genome_path, cur_sample_name, cur_tmp_dir, partition_index, open_REPCLASS):
    if os.path.exists(cur_consensus_path):
        # Step1: Running WebTE HOMOLOGY module
        file_removed = []
        #print("Thread idx:%d, Start Running WebTE HOMOLOGY module." %partition_index)
        starttime = time.time()
        #program_path = os.getcwd() + '/third-party/RepeatClassifier-2.0.1/RepeatClassifier'
        # 分类之前先在文件的最后加一条测试序列，分类后再删掉
        with open(cur_consensus_path, 'a') as f_save:
            f_save.write('>test\nTTTTTTTTTTTTTTTTTTTTTT\n')
        f_save.close()
        homology_command = 'cd '+ cur_tmp_dir + ' && RepeatClassifier -pa 1 -consensi ' + cur_consensus_path
        #print('homology_command:%s' % homology_command)
        #os.system(homology_command)
        os.system(homology_command + ' > /dev/null 2>&1')
        endtime = time.time()
        dtime = endtime - starttime
        #print("Thread idx:%d, Finish WebTE HOMOLOGY module, running time: %.4s s" % (partition_index, dtime))

        # Step2: get unclassified consensus sequence after HOMOLOGY module
        #print("get unclassified consensus sequence after HOMOLOGY module...")
        classified_consensus_path = cur_consensus_path + '.classified'
        file_removed.append(classified_consensus_path)

        if not open_REPCLASS:
            final_classified_path = cur_tmp_dir + '/consensus.fasta.final.classified'
            # os.system('mv ' + classified_consensus_path + ' ' + final_classified_path)
            # 去掉最后一行test
            names, contigs = read_fasta(classified_consensus_path)
            with open(final_classified_path, 'w') as f_save:
                for name in names:
                    if name.startswith('test'):
                        continue
                    else:
                        f_save.write('>'+name+'\n'+contigs[name]+'\n')
            f_save.close()
        else:
            classified_consensus_contignames, classified_consensus_contigs = read_fasta(classified_consensus_path)

            homology_classified_consensus_path = cur_consensus_path + '.homology.classified'
            unknown_consensus_path = cur_consensus_path + '.homology.unknown'
            file_removed.append(homology_classified_consensus_path)
            file_removed.append(unknown_consensus_path)

            known_consensus = []
            unknown_consensus = []
            for name in classified_consensus_contignames:
                name_parts = name.split('#')
                seq = classified_consensus_contigs[name]
                if name_parts[1] == 'Unknown':
                    new_name = name_parts[0]
                    unknown_consensus.append((new_name, seq))
                else:
                    known_consensus.append((name, seq))

            with open(homology_classified_consensus_path, 'w') as f_save:
                for item in known_consensus:
                    f_save.write('>' + item[0] + '\n' + item[1] + '\n')
            f_save.close()

            with open(unknown_consensus_path, 'w') as f_save:
                for item in unknown_consensus:
                    f_save.write('>' + item[0] + '\n' + item[1] + '\n')
            f_save.close()

            # Step3: Running WebTE Structure and TSD module
            print("Thread idx:%d, Start Running WebTE Structure and TSD module." %partition_index)
            starttime = time.time()
            cur_output = cur_tmp_dir + '/REPCLASS/output'
            if not os.path.exists(cur_output):
                os.makedirs(cur_output)
            #cur_user_config_dir = cur_tmp_dir + '/REPCLASS/conf'
            #set_REPCLASS_conf(cur_user_config_dir, cur_sample_name, genome_path, unknown_consensus_path, cur_output, partition_index)
            REPCLASS_home = os.getcwd() + '/third-party/REPCLASS-master/1.0.1/bin'
            user_config_path = os.getcwd() + "/third-party/REPCLASS-master/1.0.1/Myconf/userconfigure.conf"
            (GENOME_LOCATION, GENOME_FILE) = os.path.split(genome_path)
            structure_tsd_command = 'perl ' + REPCLASS_home + '/rc.pl ' + user_config_path + ' ' + cur_sample_name + \
                                    ' ' + cur_output + ' ' + GENOME_LOCATION + ' ' + GENOME_FILE + ' ' + \
                                    unknown_consensus_path + ' ' + str(partition_index)
            print('structure_tsd_command:%s' % structure_tsd_command)
            os.system(structure_tsd_command)
            endtime = time.time()
            dtime = endtime - starttime
            print("Thread idx:%d, Finish Running WebTE Structure and TSD module, running time: %.4s s" % (partition_index, dtime))

            # Step4: merge Homology classified and Structure classified
            print("Thread idx:%d, Merge Homology classified and Structure classified." % partition_index)
            REPCLASS_result_path = cur_tmp_dir + '/REPCLASS/output/Final_parsed.txt'
            marked_header = getREPCLASS_classified(REPCLASS_result_path)
            homology_classified_path = cur_tmp_dir + '/consensus.fasta.homology.classified'
            homology_unknown_path = cur_tmp_dir + '/consensus.fasta.homology.unknown'
            final_classified_path = cur_tmp_dir + '/consensus.fasta.final.classified'
            merge_classified(marked_header, homology_classified_path, homology_unknown_path, final_classified_path)


def merge_fasta(fasta_path, merged_file):
    if os.path.exists(fasta_path):
        contignames, contigs = read_fasta(fasta_path)
        with open(merged_file, 'a') as f_save:
            for name in contignames:
                seq = contigs[name]
                f_save.write('>'+name+'\n'+seq+'\n')
        f_save.close()


# if __name__ == '__main__':
#     cur_consensus_path = '/public/home/hpc194701009/repeat_detect_tools/LongRepMarker-maven/cut/TEClassTmpOutput.2021-11-23.15-46-10/1/consensus.fasta'
#     genome_path = '/public/home/hpc194701009/repeat_detect_tools/WebTE/classification/demo/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa'
#     cur_sample_name = 'dmel'
#     cur_tmp_dir = '/public/home/hpc194701009/repeat_detect_tools/LongRepMarker-maven/cut/TEClassTmpOutput.2021-11-23.15-46-10/1/'
#     partition_index = 1
#     run_classification(cur_consensus_path, genome_path, cur_sample_name, cur_tmp_dir, partition_index)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run TE classification ...')
    parser.add_argument('--sample_name', metavar='Sample Name',
                        help='input sample_name')
    parser.add_argument('--consensus', metavar='TE consensus path',
                        help='input TE consensus path')
    parser.add_argument('--genome', metavar='Genome path',
                        help='input genome path')
    parser.add_argument('--thread_num', metavar='Thread number',
                        help='input thread number')
    parser.add_argument('--split_num', metavar='Split number',
                        help='input split number')
    parser.add_argument('-o', metavar='Output Dir',
                        help='input output dir')
    args = parser.parse_args()

    sample_name = args.sample_name
    consensus_path = args.consensus
    genome_path = args.genome
    thread_num = args.thread_num
    split_num = args.split_num
    output_dir = args.o

    if not os.path.isabs(consensus_path):
        consensus_path = os.path.abspath(consensus_path)
    if genome_path != 'module' and not os.path.isabs(genome_path):
        genome_path = os.path.abspath(genome_path)
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    open_REPCLASS = False


    # Step 1: split consensus file by PET algorithm
    # partitions num equal to thread num
    partitions_num = int(split_num)
    consensus_contignames, consensus_contigs = read_fasta(consensus_path)
    data_partitions = PET(consensus_contigs.items(), partitions_num)

    # # Step 2: use thread pool to execute task parallelization
    # # delete ab-blast db file
    # db_path = os.getcwd() + '/third-party/REPCLASS-master/dependency/ab-blast/ab-blast-20200317-linux-x64/db'
    # os.system('rm -f ' + db_path + '/*')

    # create temp directory
    i = datetime.datetime.now()
    tmp_output_dir = output_dir + '/TEClassTmpOutput.' + str(i.date()) + '.' + str(i.hour) + '-' + str(i.minute) + '-' + str(i.second)
    #tmp_output_dir = output_dir + '/TEClassTmpOutput.2022-07-13.19-19-9'
    pool = multiprocessing.Pool(processes=int(thread_num))
    for partition_index, data_partition in enumerate(data_partitions):
        if len(data_partition) <= 0:
            continue
        cur_tmp_dir = tmp_output_dir + '/' + str(partition_index)
        if not os.path.exists(cur_tmp_dir):
            os.makedirs(cur_tmp_dir)
        cur_consensus_path = cur_tmp_dir + '/consensus.fasta'
        store2file(data_partition, cur_consensus_path)
        cur_sample_name = sample_name
        # print('current run %d-th species run minimap2, reference path:%s' % (species_index + 1, reference))
        pool.apply_async(run_classification, (cur_consensus_path, genome_path, cur_sample_name, cur_tmp_dir, partition_index, open_REPCLASS, ))
    pool.close()
    pool.join()

    # Step 3: merge final classified of each thread
    final_classified_path = consensus_path + '.classified'
    final_tmpBlastX_path = consensus_path + '.tmpBlastXResults.out.bxsummary'
    final_tmpBlastn_path = consensus_path + '.tmpBlastnResults.out'
    if os.path.exists(final_classified_path):
        os.system('rm -f '+ final_classified_path)
    if os.path.exists(final_tmpBlastX_path):
        os.system('rm -f '+ final_tmpBlastX_path)
    if os.path.exists(final_tmpBlastn_path):
        os.system('rm -f '+ final_tmpBlastn_path)
    for partition_index in range(partitions_num):
        cur_tmp_dir = tmp_output_dir + '/' + str(partition_index)
        cur_classified_path = cur_tmp_dir + '/consensus.fasta.final.classified'
        if os.path.exists(cur_classified_path):
            merge_command = 'cat ' + cur_classified_path + ' >> ' + final_classified_path
            #print(merge_command)
            os.system(merge_command)

        cur_tmpBlastxOutput = cur_tmp_dir + '/tmpBlastXResults.out.bxsummary'
        if os.path.exists(cur_tmpBlastxOutput):
            merge_command = 'cat ' + cur_tmpBlastxOutput + ' >> ' + final_tmpBlastX_path
            #print(merge_command)
            os.system(merge_command)

        cur_tmpBlastnOutput = cur_tmp_dir + '/blastn.out'
        if os.path.exists(cur_tmpBlastnOutput):
            merge_command = 'cat ' + cur_tmpBlastnOutput + ' >> ' + final_tmpBlastn_path
            #print(merge_command)
            os.system(merge_command)
            #merge_fasta(cur_classified_path, final_classified_path)





