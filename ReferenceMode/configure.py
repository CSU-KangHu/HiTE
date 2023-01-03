import os
import json
import sys

def valid_path(path):
    if path == '' or path is None or not os.path.exists(path):
        return False
    else:
        return True


def modify_RepModel(RepModelConfig_path, repeatmasker_home, rmblast_home, trf_program):
    RepModelConfig_lines = []

    modifing_repeatmasker = False
    modifing_rmblast = False
    modifing_trf = False
    modified_repeatmasker = False
    modified_rmblast = False
    modified_trf = False

    with open(RepModelConfig_path, 'r') as f_r:
        for line in f_r.readlines():
            if not modified_repeatmasker:
                if not modifing_repeatmasker:
                    if line.split('=>')[0].strip() == "'REPEATMASKER_DIR'":
                        modifing_repeatmasker = True
                else:
                    if line.split('=>')[0].strip() == "'value'":
                        # modify
                        new_line = line.split('=>')[0] + "=> '" + repeatmasker_home + "'\n"
                        modified_repeatmasker = True
                        line = new_line

            if not modified_rmblast:
                if not modifing_rmblast:
                    if line.split('=>')[0].strip() == "'RMBLAST_DIR'":
                        modifing_rmblast = True
                else:
                    if line.split('=>')[0].strip() == "'value'":
                        # modify
                        new_line = line.split('=>')[0] + "=> '" + rmblast_home + "/bin'\n"
                        modified_rmblast = True
                        line = new_line

            if not modified_trf:
                if not modifing_trf:
                    if line.split('=>')[0].strip() == "'TRF_PRGM'":
                        modifing_trf = True
                else:
                    if line.split('=>')[0].strip() == "'value'":
                        # modify
                        new_line = line.split('=>')[0] + "=> '" + trf_program + "'\n"
                        modified_trf = True
                        line = new_line

            RepModelConfig_lines.append(line)

    with open(RepModelConfig_path, 'w') as f_w:
        f_w.writelines(RepModelConfig_lines)


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

    with open(repclass_conf_path, 'w') as f_save:
        f_save.write(new_str)

if __name__ == '__main__':
    param_config_path = os.getcwd() + "/ParamConfig.json"
    # read param config
    with open(param_config_path, 'r') as load_f:
        param = json.load(load_f)
    load_f.close()

    # 1. get RepeatClassifier path
    repeatclassifier_path = os.getcwd() + "/classification/third-party/RepeatClassifier-2.0.1/RepeatClassifier"
    if not os.path.exists(repeatclassifier_path):
        print('error:')
        print('----can not find RepeatClassifier in %s' %repeatclassifier_path)
        sys.exit(-1)

    repeatmasker_home = param['RepeatMasker_Home']
    if not os.path.exists(repeatmasker_home+'/RepeatMasker'):
        print('error:')
        print('----can not find RepeatMasker in %s directory' %repeatmasker_home)
        sys.exit(-1)

    rmblast_home = param['RMBlast_Home']
    if not os.path.exists(rmblast_home+'/bin/makeblastdb'):
        print('error:')
        print('----can not find makeblastdb in %s directory' %rmblast_home)
        sys.exit(-1)

    genome_tools_home = param['Genome_Tools_Home']
    if not os.path.exists(genome_tools_home+'/bin/gt'):
        print('error:')
        print('----can not find gt in %s directory' %genome_tools_home)
        sys.exit(-1)

    ltr_retriever_home = param['LTR_retriever_Home']
    if not os.path.exists(ltr_retriever_home+'/LTR_retriever'):
        print('error:')
        print('----can not find LTR_retriever in %s directory' %ltr_retriever_home)
        sys.exit(-1)

    LTR_finder_parallel_home = param['LTR_finder_parallel_Home'] + "/LTR_FINDER_parallel"
    if not os.path.exists(LTR_finder_parallel_home):
        print('error:')
        print('----can not find LTR_finder_parallel in %s' % LTR_finder_parallel_home)
        sys.exit(-1)

    trf_program = os.getcwd() + "/tools/trf409.linux64"
    if not os.path.exists(trf_program):
        print('error:')
        print('----can not find trf in %s' %trf_program)
        sys.exit(-1)

    EAHelitron_program = param['EAHelitron'] + "/EAHelitron"
    if not os.path.exists(EAHelitron_program):
        print('error:')
        print('----can not find EAHelitron in %s' % EAHelitron_program)
        sys.exit(-1)


    # 2. modify RepModelConfig.pm in RepeatClassifier directory
    RepModelConfig_path = os.getcwd() + "/classification/third-party/RepeatClassifier-2.0.1/RepModelConfig.pm"
    modify_RepModel(RepModelConfig_path, repeatmasker_home, rmblast_home, trf_program)

    # 3. add executable
    os.system('chmod +x ' + repeatclassifier_path)

    #print('Start set REPCLASS...')
    # 1. set repclass.conf
    #repclass_conf = {}
    #repclass_conf_path = os.getcwd() + "/third-party/REPCLASS-master/1.0.1/conf/repclass.conf"
    #REPCLASS = os.getcwd() + "/third-party/REPCLASS-master/1.0.1"
    #PROGRAMS = os.getcwd() + "/third-party/REPCLASS-master/dependency/EMBOSS-6.5.7/build/bin"
    #WUBLAST = os.getcwd() + "/third-party/REPCLASS-master/dependency/ab-blast/ab-blast-20200317-linux-x64"
    #repclass_conf['REPCLASS'] = REPCLASS
    #repclass_conf['PROGRAMS'] = PROGRAMS
    #repclass_conf['WUBLAST'] = WUBLAST
    #modify_REPCLASS_conf(repclass_conf_path, repclass_conf)

    # 1. set repclass.conf
    #repclass_conf = {}
    #repclass_conf_path = os.getcwd() + "/third-party/REPCLASS-master/1.0.1/conf/repclass.conf"
    #REPCLASS = os.getcwd() + "/third-party/REPCLASS-master/1.0.1"
    #PROGRAMS = os.getcwd() + "/third-party/REPCLASS-master/dependency/EMBOSS-6.5.7/build/bin"
    #WUBLAST = os.getcwd() + "/third-party/REPCLASS-master/dependency/ab-blast/ab-blast-20200317-linux-x64"
    #repclass_conf['REPCLASS'] = REPCLASS
    #repclass_conf['PROGRAMS'] = PROGRAMS
    #repclass_conf['WUBLAST'] = WUBLAST
    #modify_REPCLASS_conf(repclass_conf_path, repclass_conf)

    # 2. set userconfigure.conf
    #userconfigure_conf = {}
    #userconfigure_conf_path = os.getcwd() + '/third-party/REPCLASS-master/1.0.1/Myconf/userconfigure.conf'
    #repclass_conf_path = os.getcwd() + "/third-party/REPCLASS-master/1.0.1/conf/repclass.conf"
    #TEMPCHOICE = os.getcwd() + '/third-party/REPCLASS-master/1.0.1/Myconf/mytempchoice'
    #REPCLASS_CONF = repclass_conf_path

    #userconfigure_conf['TEMPCHOICE'] = TEMPCHOICE
    #userconfigure_conf['REPCLASS_CONF'] = REPCLASS_CONF
    #modify_REPCLASS_conf(userconfigure_conf_path, userconfigure_conf)

    #print('Finish set REPCLASS...')

    print("Congratulation! all configuration pass.")
