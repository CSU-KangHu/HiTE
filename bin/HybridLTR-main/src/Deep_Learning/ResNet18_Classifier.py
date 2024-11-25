#-- coding: UTF-8 --
import argparse
import json
import os
import sys
import time

import numpy as np
from matplotlib import pyplot as plt



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
from ResNet18 import ResNet18
import tensorflow as tf
from PIL import Image
from Util import find_files_recursively

class Classifier:
    def __init__(self, model_path):
        # set input image parameters
        image_size = (100, 200)
        channels = 3
        num_classes = 2

        self.model = ResNet18(num_classes)
        self.model.build(input_shape=(None, image_size[1], image_size[0], channels))
        self.model.load_weights(model_path)


    def predict(self, test_dir, y, row_nums, matrix_files):
        image_shape = config.image_shape
        image_size = (image_shape[0], image_shape[1])
        channels = 3
        num_classes = 2
        batch_size = 1

        seq_names = []
        X = []
        # 遍历目录中的所有图片文件
        for filename in os.listdir(test_dir):
            seq_names.append(filename.split('.')[0])
            if filename.endswith(".jpg") or filename.endswith(".png"):  # 检查文件扩展名
                file_path = os.path.join(test_dir, filename)
                image = Image.open(file_path)
                # plt.imshow(image)
                # plt.show()
                image = image.resize(image_size)
                image_array = np.array(image) / 255.0  # 归一化到0-1之间
                X.append(image_array)
        X = np.array(X)


        y_pred = self.model.predict(X)
        print(y_pred)
        y_pred = np.argmax(np.round(y_pred), axis=1)
        print(y_pred)

        output_path = config.work_dir + '/is_LTR_deep.txt'
        with open(output_path, 'w') as f_save:
            for i in range(len(y_pred)):
                ltr_name = seq_names[i]
                f_save.write(ltr_name + '\t' + str(y_pred[i]) + '\n')


def get_image_from_matrix(matrix_file):
    image_shape = config.image_shape
    row_num = image_shape[0]
    col_num = image_shape[1]
    lines = []
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').replace('\t', '')
            lines.append(line)
    if len(lines) < 4:
        return None

    matrix = [['-'] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

    image = np.zeros((row_num, col_num, image_shape[2]))
    base_value = {'-': (255, 255, 255), 'A': (59, 197, 53), 'T': (245, 97, 97), 'C': (88, 221, 255),
                  'G': (185, 80, 250)}
    for col_index in range(col_num):
        for row_index in range(row_num):
            base = matrix[row_index][col_index]
            if base not in base_value:
                cur_val = (255, 255, 255)
            else:
                cur_val = base_value[base]
            image[row_index, col_index] = cur_val
    return image
def convert_matrix2jpg(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    positive_dir = input_dir + '/positive'
    positive_output_dir = output_dir + '/positive'
    if not os.path.exists(positive_output_dir):
        os.makedirs(positive_output_dir)

    file_extension = '.matrix'
    all_positive_matrix_files = find_files_recursively(positive_dir, file_extension)
    for matrix_file in all_positive_matrix_files:
        image = get_image_from_matrix(matrix_file)
        if image is not None:
            name = os.path.basename(matrix_file)
            image_path = positive_output_dir + '/' + name + '.jpg'
            image = Image.fromarray(np.uint8(image))
            image.save(image_path)

    negative_dir = input_dir + '/negative'
    negative_output_dir = output_dir + '/negative'
    if not os.path.exists(negative_output_dir):
        os.makedirs(negative_output_dir)
    all_negative_matrix_files = find_files_recursively(negative_dir, file_extension)
    for matrix_file in all_negative_matrix_files:
        image = get_image_from_matrix(matrix_file)
        if image is not None:
            name = os.path.basename(matrix_file)
            image_path = negative_output_dir + '/' + name + '.jpg'
            image = Image.fromarray(np.uint8(image))
            image.save(image_path)

def main():
    showToolName()

    # 1.parse args
    describe_info = '########################## NeuralTE, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data_dir', required=True, metavar='data_dir', help='Input directory contains positive and negative data')
    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='Output directory, store temporary files')
    parser.add_argument('--model_path', required=True, metavar='model_path', help='Input the path of trained model, absolute path.')

    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')



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

    # reload GPU config
    get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

    starttime1 = time.time()
    # Instantiate the DataProcessor class
    data_processor = DataProcessor()
    # load data
    test_dir = config.work_dir + '/negative'
    if os.path.exists(test_dir):
        os.system('rm -rf ' + test_dir)
    X, y, row_nums, matrix_files = data_processor.load_predict_data(data_dir)

    print(X.shape, y.shape)
    endtime1 = time.time()
    dtime1 = endtime1 - starttime1
    print("Running time of DataProcessor: %.8s s" % (dtime1))


    starttime2 = time.time()
    # Instantiate the Classifier class
    classifier = Classifier(model_path=model_path)
    # start prediction
    classifier.predict(test_dir, y, row_nums, matrix_files)
    endtime2 = time.time()
    dtime2 = endtime2 - starttime2
    print("Running time of model Classifier: %.8s s" % (dtime2))


    endtime3 = time.time()
    dtime = endtime3 - starttime1
    print("Running time of total NeuralTE-Classifier: %.8s s" % (dtime))

if __name__ == '__main__':
    main()