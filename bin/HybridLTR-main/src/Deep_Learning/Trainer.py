#-- coding: UTF-8 --
import argparse
import json
import os
import sys
import time

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "../..")
sys.path.append(configs_folder)

import numpy as np
from keras.utils import np_utils
from configs import config, gpu_config
from CNN_Model_image import CNN_Model
from LSTM_Model import LSTM_Model
from keras.callbacks import EarlyStopping
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

        # cnn_model = CNN_Model(config.X_feature_len, config.class_num)
        # model = cnn_model.build_model(cnn_num_convs, cnn_filters_array)

        # LSTM 模型
        lstm_model = LSTM_Model(config.lstm_seq_len, config.lstm_features_num, config.class_num)
        model = lstm_model.build_model(config.lstm_num_layers, config.lstm_hidden_size, config.lstm_dropout)


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
    parser.add_argument('--data_dir', required=True, metavar='data_dir', help='Input directory contains positive and negative data')
    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='Output directory, store temporary files')

    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')

    parser.add_argument('--threads', metavar='thread_num', help='Input thread num, default = [ ' + str(config.threads) + ' ]')

    parser.add_argument('--cnn_num_convs', metavar='cnn_num_convs', help='The number of CNN convolutional layers. default = [ ' + str(config.cnn_num_convs) + ' ]')
    parser.add_argument('--cnn_filters_array', metavar='cnn_filters_array', help='The number of filters in each CNN convolutional layer. default = [ ' + str(config.cnn_filters_array) + ' ]')
    parser.add_argument('--cnn_kernel_sizes_array', metavar='cnn_kernel_sizes_array', help='The kernel size in each of CNN convolutional layer. default = [ ' + str(config.cnn_kernel_sizes_array) + ' ]')
    parser.add_argument('--cnn_dropout', metavar='cnn_dropout', help='The threshold of CNN Dropout. default = [ ' + str(config.cnn_dropout) + ' ]')
    parser.add_argument('--batch_size', metavar='batch_size', help='The batch size in training model. default = [ ' + str(config.batch_size) + ' ]')
    parser.add_argument('--epochs', metavar='epochs', help='The number of epochs in training model. default = [ ' + str(config.epochs) + ' ]')
    parser.add_argument('--use_checkpoint', metavar='use_checkpoint',  help='Whether to use breakpoint training. 1: true, 0: false. The model will continue training from the last failed parameters to avoid training from head. default = [ ' + str(config.use_checkpoint) + ' ]')


    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    start_gpu_num = args.start_gpu_num
    use_gpu_num = args.use_gpu_num

    threads = args.threads
    cnn_num_convs = args.cnn_num_convs
    cnn_filters_array = args.cnn_filters_array
    cnn_kernel_sizes_array = args.cnn_kernel_sizes_array
    cnn_dropout = args.cnn_dropout
    batch_size = args.batch_size
    epochs = args.epochs
    use_checkpoint = args.use_checkpoint

    if out_dir is not None:
        config.work_dir = out_dir

    if start_gpu_num is not None:
        gpu_config.start_gpu_num = int(start_gpu_num)
    if use_gpu_num is not None:
        gpu_config.use_gpu_num = int(use_gpu_num)


    if threads is not None:
        config.threads = int(threads)

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


    params = {}
    params['data_dir'] = data_dir
    params['out_dir'] = out_dir
    showTrainParams(params)

    positive_dir = data_dir + '/positive'
    negative_dir = data_dir + '/negative'

    # re-compute feature length
    row_num = 100
    col_num = 100
    X_feature_len = get_feature_len(row_num, col_num)
    config.X_feature_len = X_feature_len

    config.lstm_seq_len = col_num + 1
    config.lstm_features_num = 5

    # reload GPU config
    get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

    starttime1 = time.time()
    # Instantiate the DataProcessor class
    data_processor = DataProcessor()
    # load data
    X, y, row_nums, matrix_files = data_processor.load_data(positive_dir, negative_dir)
    print(X.shape, y.shape)
    endtime1 = time.time()
    dtime1 = endtime1 - starttime1
    print("Running time of DataProcessor: %.8s s" % (dtime1))

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