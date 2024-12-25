#-- coding: UTF-8 --
import argparse
import json
import os
import sys
import time

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

import numpy as np
from keras.utils import np_utils
from configs import config, gpu_config
from CNN_Model import CNN_Model
from DataProcessor import DataProcessor
from utils.show_util import showToolName, showTrainParams
from utils.data_util import get_feature_len, get_gpu_config
import datetime



class Trainer:
    def __int__(self):
        pass

    def train(self, X, y, cnn_num_convs, cnn_filters_array):
        y_one_hot = np_utils.to_categorical(y, int(config.class_num))
        y_one_hot = np.array(y_one_hot)

        cnn_model = CNN_Model(config.X_feature_len, config.class_num)
        model = cnn_model.build_model(cnn_num_convs, cnn_filters_array)

        # model training
        model.fit(X, y_one_hot, batch_size=config.batch_size, epochs=config.epochs, verbose=1)
        # save model
        i = datetime.datetime.now()
        time_str = str(i.date()) + '.' + str(i.hour) + '-' + str(i.minute) + '-' + str(i.second)
        model_path = config.project_dir + '/models/' + f'model_' + time_str + '.h5'
        model.save(model_path)
        return model_path


def main():
    showToolName()

    # 1.parse args
    describe_info = '########################## NeuralTE, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data', required=True, metavar='data', help='Input fasta file used to train model, header format: seq_name\tlabel\tspecies_name, refer to "data/train.example.fa" for example.')
    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='Output directory, store temporary files')
    parser.add_argument('--use_TSD', required=True, metavar='use_TSD', help='Whether to use TSD features, 1: true, 0: false. default = [ ' + str(config.use_TSD) + ' ]')
    parser.add_argument('--is_train', required=True, metavar='is_train', help='Enable train mode, 1: true, 0: false. default = [ ' + str(config.is_train) + ' ]')
    parser.add_argument('--is_predict', required=True, metavar='is_predict', help='Enable prediction mode, 1: true, 0: false. default = [ ' + str(config.is_predict) + ' ]')

    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')
    parser.add_argument('--only_preprocess', metavar='only_preprocess', help='Whether to only perform data preprocessing, 1: true, 0: false.')
    parser.add_argument('--keep_raw', metavar='keep_raw', help='Whether to retain the raw input sequence, 1: true, 0: false; only save species having TSDs. default = [ ' + str(config.keep_raw) + ' ]')
    parser.add_argument('--genome', metavar='genome', help='Genome path, use to search for TSDs')
    parser.add_argument('--use_kmers', metavar='use_kmers', help='Whether to use kmers features, 1: true, 0: false. default = [ ' + str(config.use_kmers) + ' ]')
    parser.add_argument('--use_terminal', metavar='use_terminal', help='Whether to use LTR, TIR terminal features, 1: true, 0: false. default = [ ' + str(config.use_terminal) + ' ]')
    parser.add_argument('--use_minority', metavar='use_minority', help='Whether to use minority features, 1: true, 0: false. default = [ ' + str(config.use_minority) + ' ]')
    parser.add_argument('--use_domain', metavar='use_domain', help='Whether to use domain features, 1: true, 0: false. default = [ ' + str(config.use_domain) + ' ]')
    parser.add_argument('--use_ends', metavar='use_ends', help='Whether to use 5-bp terminal ends features, 1: true, 0: false. default = [ ' + str(config.use_ends) + ' ]')
    parser.add_argument('--threads', metavar='thread_num', help='Input thread num, default = [ ' + str(config.threads) + ' ]')
    parser.add_argument('--internal_kmer_sizes', metavar='internal_kmer_sizes', help='The k-mer size used to convert internal sequences to k-mer frequency features, default = [ ' + str(config.internal_kmer_sizes) + ' MB ]')
    parser.add_argument('--terminal_kmer_sizes', metavar='terminal_kmer_sizes', help='The k-mer size used to convert terminal sequences to k-mer frequency features, default = [ ' + str(config.terminal_kmer_sizes) + ' ]')
    parser.add_argument('--cnn_num_convs', metavar='cnn_num_convs', help='The number of CNN convolutional layers. default = [ ' + str(config.cnn_num_convs) + ' ]')
    parser.add_argument('--cnn_filters_array', metavar='cnn_filters_array', help='The number of filters in each CNN convolutional layer. default = [ ' + str(config.cnn_filters_array) + ' ]')
    parser.add_argument('--cnn_kernel_sizes_array', metavar='cnn_kernel_sizes_array', help='The kernel size in each of CNN convolutional layer. default = [ ' + str(config.cnn_kernel_sizes_array) + ' ]')
    parser.add_argument('--cnn_dropout', metavar='cnn_dropout', help='The threshold of CNN Dropout. default = [ ' + str(config.cnn_dropout) + ' ]')
    parser.add_argument('--batch_size', metavar='batch_size', help='The batch size in training model. default = [ ' + str(config.batch_size) + ' ]')
    parser.add_argument('--epochs', metavar='epochs', help='The number of epochs in training model. default = [ ' + str(config.epochs) + ' ]')
    parser.add_argument('--use_checkpoint', metavar='use_checkpoint',  help='Whether to use breakpoint training. 1: true, 0: false. The model will continue training from the last failed parameters to avoid training from head. default = [ ' + str(config.use_checkpoint) + ' ]')


    args = parser.parse_args()

    data_path = args.data
    out_dir = args.out_dir
    genome = args.genome
    use_kmers = args.use_kmers
    use_terminal = args.use_terminal
    use_TSD = args.use_TSD
    use_minority = args.use_minority
    use_domain = args.use_domain
    use_ends = args.use_ends
    is_train = args.is_train
    only_preprocess = args.only_preprocess

    start_gpu_num = args.start_gpu_num
    use_gpu_num = args.use_gpu_num
    keep_raw = args.keep_raw
    is_predict = args.is_predict
    threads = args.threads
    internal_kmer_sizes = args.internal_kmer_sizes
    terminal_kmer_sizes = args.terminal_kmer_sizes
    cnn_num_convs = args.cnn_num_convs
    cnn_filters_array = args.cnn_filters_array
    cnn_kernel_sizes_array = args.cnn_kernel_sizes_array
    cnn_dropout = args.cnn_dropout
    batch_size = args.batch_size
    epochs = args.epochs
    use_checkpoint = args.use_checkpoint

    if out_dir is not None:
        config.work_dir = out_dir
    if use_kmers is not None:
        config.use_kmers = int(use_kmers)
    if use_terminal is not None:
        config.use_terminal = int(use_terminal)
    if use_TSD is not None:
        config.use_TSD = int(use_TSD)
    if use_minority is not None:
        config.use_minority = int(use_minority)
    if use_domain is not None:
        config.use_domain = int(use_domain)
    if use_ends is not None:
        config.use_ends = int(use_ends)
    if is_train is not None:
        config.is_train = int(is_train)

    if start_gpu_num is not None:
        gpu_config.start_gpu_num = int(start_gpu_num)
    if use_gpu_num is not None:
        gpu_config.use_gpu_num = int(use_gpu_num)
    if keep_raw is not None:
        config.keep_raw = int(keep_raw)
    if only_preprocess is not None:
        config.only_preprocess = int(only_preprocess)
    if is_predict is not None:
        config.is_predict = int(is_predict)
    if threads is not None:
        config.threads = int(threads)
    if internal_kmer_sizes is not None:
        config.internal_kmer_sizes = json.loads(internal_kmer_sizes)
    if terminal_kmer_sizes is not None:
        config.terminal_kmer_sizes = json.loads(terminal_kmer_sizes)
    if cnn_num_convs is not None:
        config.cnn_num_convs = int(cnn_num_convs)
    if cnn_filters_array is not None:
        config.cnn_filters_array = json.loads(cnn_filters_array)
    if cnn_kernel_sizes_array is not None:
        config.cnn_kernel_sizes_array = json.loads(cnn_kernel_sizes_array)
    if cnn_dropout is not None:
        config.cnn_dropout = float(cnn_dropout)
    if batch_size is not None:
        config.batch_size = int(batch_size)
    if epochs is not None:
        config.epochs = int(epochs)
    if use_checkpoint is not None:
        config.use_checkpoint = int(use_checkpoint)

    if genome is not None:
        os.makedirs(config.work_dir, exist_ok=True)
        genome_info_path = config.work_dir + '/genome.info'
        if str(genome).__contains__('genome.info'):
            os.system('cp ' + genome + ' ' + genome_info_path)
        else:
            with open(genome_info_path, 'w') as f_save:
                f_save.write('#Scientific Name\tGenome Path\tIs Plant\n')
                f_save.write('Unknown\t'+genome+'\t'+str(config.is_plant)+'\n')

    params = {}
    params['data_path'] = data_path
    params['out_dir'] = out_dir
    params['genome'] = genome
    showTrainParams(params)

    # re-compute feature length
    X_feature_len = get_feature_len()
    config.X_feature_len = X_feature_len

    # reload GPU config
    get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

    starttime1 = time.time()
    # Instantiate the DataProcessor class
    data_processor = DataProcessor()
    # load data
    X, y, seq_names, data_path = data_processor.load_data(config.internal_kmer_sizes, config.terminal_kmer_sizes, data_path)
    print(X.shape, y.shape)
    endtime1 = time.time()
    dtime1 = endtime1 - starttime1
    print("Running time of DataProcessor: %.8s s" % (dtime1))

    if not config.only_preprocess:
        starttime2 = time.time()
        # Instantiate the Trainer class
        trainer = Trainer()
        # start training
        model_path = trainer.train(X, y, config.cnn_num_convs, config.cnn_filters_array)
        print('Trained model is stored in:', model_path)
        endtime2 = time.time()
        dtime2 = endtime2 - starttime2
        print("Running time of model Trainer: %.8s s" % (dtime2))

    endtime3 = time.time()
    dtime = endtime3 - starttime1
    print("Running time of total NeuralTE-Trainer: %.8s s" % (dtime))

if __name__ == '__main__':
    main()