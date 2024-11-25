#-- coding: UTF-8 --
import argparse
import json
import os
import sys
import time

import numpy as np

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "../..")
sys.path.append(configs_folder)

from utils.show_util import showToolName, showTestParams
from keras.models import load_model
from configs import config
from configs import gpu_config
from DataProcessor import DataProcessor
from utils.evaluate_util import get_metrics, correct_using_minority, get_metrics_by_label
from utils.data_util import get_feature_len, get_gpu_config
from sklearn.preprocessing import StandardScaler

class Classifier:
    def __init__(self, model_path):
        self.model = load_model(model_path)
        #self.model = load_model(model_path, compile=False)

    def predict(self, X, y, row_nums, matrix_files):
        # 预测概率
        y_pred = self.model.predict(X)
        y_pred = np.argmax(np.round(y_pred), axis=1)

        output_path = config.work_dir + '/is_LTR_deep.txt'
        with open(output_path, 'w') as f_save:
            for i in range(len(y_pred)):
                ltr_name = os.path.basename(matrix_files[i]).replace('.matrix', '')
                f_save.write(ltr_name + '\t' + str(y_pred[i]) + '\n')

    def predict_hybrid(self, X1, X2, y, row_nums, matrix_files):

        # 预测概率
        y_pred = self.model.predict([X1, X2])
        y_pred = np.argmax(np.round(y_pred), axis=1)

        output_path = config.work_dir + '/is_LTR_deep.txt'
        with open(output_path, 'w') as f_save:
            for i in range(len(y_pred)):
                ltr_name = os.path.basename(matrix_files[i]).replace('.matrix', '')
                f_save.write(ltr_name + '\t' + str(y_pred[i]) + '\n')


def main():
    showToolName()

    # 1.parse args
    describe_info = '########################## NeuralTE, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data_dir', required=True, metavar='data_dir', help='Input directory contains positive and negative data')
    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='Output directory, store temporary files')

    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')

    parser.add_argument('--model_path', metavar='model_path', help='Input the path of trained model, absolute path.')

    parser.add_argument('--threads', metavar='thread_num', help='Input thread num, default = [ ' + str(config.threads) + ' ]')

    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    model_path = args.model_path

    start_gpu_num = args.start_gpu_num
    use_gpu_num = args.use_gpu_num

    threads = args.threads


    if out_dir is not None:
        config.work_dir = out_dir

    if start_gpu_num is not None:
        gpu_config.start_gpu_num = int(start_gpu_num)
    if use_gpu_num is not None:
        gpu_config.use_gpu_num = int(use_gpu_num)
    if threads is not None:
        config.threads = int(threads)

    params = {}
    params['data_dir'] = data_dir
    params['out_dir'] = out_dir
    params['model_path'] = model_path
    showTestParams(params)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # re-compute feature length
    row_num = 100
    col_num = 100
    # X_feature_len = 0
    # for kmer_size in config.kmer_sizes:
    #     X_feature_len += 100 - kmer_size + 1
    X_feature_len = 0
    for kmer_size in config.kmer_sizes:
        X_feature_len += pow(5, kmer_size) - 1
    # X_feature_len += 5 * 100
    config.X_feature_len = X_feature_len

    config.lstm_seq_len = 100
    config.lstm_features_num = 5

    # reload GPU config
    get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

    starttime1 = time.time()
    # Instantiate the DataProcessor class
    data_processor = DataProcessor()
    # # load data
    # X, y, row_nums, matrix_files = data_processor.load_predict_data_hybrid(data_dir)
    # print(X.shape,y.shape)

    # Load the data
    X1, X2, y, row_nums, matrix_files = data_processor.load_predict_data_hybrid(data_dir)
    print(X1.shape, X2.shape, y.shape)
    endtime1 = time.time()
    dtime1 = endtime1 - starttime1
    print("Running time of DataProcessor: %.8s s" % (dtime1))

    starttime2 = time.time()
    # Instantiate the Classifier class
    classifier = Classifier(model_path=model_path)
    # start prediction
    # classifier.predict(X, y, row_nums, matrix_files)

    classifier.predict_hybrid(X1, X2, y, row_nums, matrix_files)
    endtime2 = time.time()
    dtime2 = endtime2 - starttime2
    print("Running time of model Classifier: %.8s s" % (dtime2))


    endtime3 = time.time()
    dtime = endtime3 - starttime1
    print("Running time of total NeuralTE-Classifier: %.8s s" % (dtime))

if __name__ == '__main__':
    main()