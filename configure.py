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
    f_r.close()
    with open(RepModelConfig_path, 'w') as f_w:
        f_w.writelines(RepModelConfig_lines)
    f_w.close()


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
    f_save.close()

if __name__ == '__main__':
    param_config_path = os.getcwd() + "/ParamConfig.json"
    # read param config
    with open(param_config_path, 'r') as load_f:
        param = json.load(load_f)
    load_f.close()

    found_error = False
    found_warn = False

    # 1. get RepeatClassifier path
    genome_tools_home = param['Genome_Tools_Home']
    if not os.path.exists(genome_tools_home+'/bin/gt'):
        print('error:')
        print('----can not find gt in %s directory.' %genome_tools_home)
        found_error = True

    ltr_retriever_home = param['LTR_retriever_Home']
    if not os.path.exists(ltr_retriever_home+'/LTR_retriever'):
        print('error:')
        print('----can not find LTR_retriever in %s directory.' %ltr_retriever_home)
        found_error = True

    rmblast_home = param['RMBlast_Home']
    if not os.path.exists(rmblast_home + '/bin/makeblastdb'):
        print('error:')
        print('----can not find makeblastdb in %s directory.' % rmblast_home)
        found_error = True

    repeatmasker_home = param['RepeatMasker_Home']
    if not os.path.exists(repeatmasker_home + '/RepeatMasker'):
        print('warning:')
        print('----can not find RepeatMasker in %s directory.' % repeatmasker_home)
        found_warn = True

    repeatModeler_Home = param['RepeatModeler_Home']
    if not os.path.exists(repeatModeler_Home + '/RepeatClassifier'):
        print('warning:')
        print('----can not find RepeatClassifier in %s directory.' % repeatModeler_Home)
        found_warn = True


    if found_error:
        print('Configuration error. Please enter the correct tool paths in ParamConfig.json.')
        sys.exit(-1)
    elif found_warn:
        print('Configuration warning. RepeatModeler or RepeatMasker is not configured correctly, please run HiTE with the \'--classified 0\' option.')

    # # 2. modify RepModelConfig.pm in RepeatClassifier directory
    # RepModelConfig_path = os.getcwd() + "/classification/third-party/RepeatClassifier-2.0.1/RepModelConfig.pm"
    # modify_RepModel(RepModelConfig_path, repeatmasker_home, rmblast_home, trf_program)

    # 3. add executable
    #os.system('cd ' + os.getcwd() + '/tools && chmod +x ./*')

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
