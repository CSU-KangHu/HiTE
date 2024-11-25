# This program loads data and determines whether data needs respective preprocessing
# based on which features are used in the configuration file.
import os
import sys

from sklearn.preprocessing import StandardScaler

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "../..")
sys.path.append(configs_folder)

from configs import config
from concurrent.futures import ProcessPoolExecutor, as_completed

from utils.data_util import read_fasta_v1, generate_TSD_info, generate_terminal_info, \
    load_repbase_with_TSD, generate_feature_mats, read_fasta, connect_LTR, store_fasta, generate_minority_info, \
    generate_predict_feature_mats, generate_feature_mats_hybrid, generate_predict_feature_mats_hybrid


class DataProcessor:
    def __init__(self):
        self.project_dir = config.project_dir
        self.tool_dir = self.project_dir + '/tools'
        self.work_dir = config.work_dir
        self.threads = config.threads

        self.ex = ProcessPoolExecutor(self.threads)


    def load_data(self, positive_dir, negative_dir):
        # Copy input files to the working directory
        if not os.path.exists(positive_dir):
            print('Positive directory not exist: ' + positive_dir)
            exit(-1)
        if not os.path.exists(negative_dir):
            print('Negative directory not exist: ' + negative_dir)
            exit(-1)

        X, Y, row_nums, matrix_files = generate_feature_mats(positive_dir, negative_dir, self.ex)

        # Reshape data into the format accepted by the model

        X = X.reshape(X.shape[0], config.X_feature_len, 1)
        X = X.astype('float64')
        # image_shape = config.image_shape
        # X = X.reshape(X.shape[0], image_shape[0], image_shape[1], image_shape[2])
        # X2 = X2.reshape(X2.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X2 = X2.astype('float64')
        # X = X.reshape(X.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X = X.astype('float64')
        return X, Y, row_nums, matrix_files

    def load_data_hybrid(self, positive_dir, negative_dir):
        # Copy input files to the working directory
        if not os.path.exists(positive_dir):
            print('Positive directory not exist: ' + positive_dir)
            exit(-1)
        if not os.path.exists(negative_dir):
            print('Negative directory not exist: ' + negative_dir)
            exit(-1)

        X1, X2, Y, row_nums, matrix_files = generate_feature_mats_hybrid(positive_dir, negative_dir, self.ex)

        scaler = StandardScaler()
        # 训练集归一化
        X1 = scaler.fit_transform(X1)
        # 验证集归一化，注意我们使用训练集的scaler来转换验证集
        X2 = scaler.transform(X2)

        # Reshape data into the format accepted by the model
        X1 = X1.reshape(X1.shape[0], config.X_feature_len, 1)
        X1 = X1.astype('float64')
        X2 = X2.reshape(X2.shape[0], config.X_feature_len, 1)
        X2 = X2.astype('float64')

        # X1 = X1.reshape(X1.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X1 = X1.astype('float64')
        # X2 = X2.reshape(X2.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X2 = X2.astype('float64')
        return X1, X2, Y, row_nums, matrix_files

    def load_predict_data(self, predict_dir):
        # Copy input files to the working directory
        if not os.path.exists(predict_dir):
            print('Predict_dir directory not exist: ' + predict_dir)
            exit(-1)

        X, Y, row_nums, matrix_files = generate_predict_feature_mats(predict_dir, self.ex)

        # Reshape data into the format accepted by the model
        # X = X.reshape(X.shape[0], config.X_feature_len, 1)
        # X = X.astype('float64')

        # image_shape = config.image_shape
        # X = X.reshape(X.shape[0], image_shape[0], image_shape[1], image_shape[2])

        # X2 = X2.reshape(X2.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X2 = X2.astype('float64')
        # X = X.reshape(X.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X = X.astype('float64')
        return X, Y, row_nums, matrix_files

    def load_predict_data_hybrid(self, predict_dir):
        # Copy input files to the working directory
        if not os.path.exists(predict_dir):
            print('Predict_dir directory not exist: ' + predict_dir)
            exit(-1)

        X1, X2, Y, row_nums, matrix_files = generate_predict_feature_mats_hybrid(predict_dir, self.ex)

        scaler = StandardScaler()
        # 训练集归一化
        X1 = scaler.fit_transform(X1)
        # 验证集归一化，注意我们使用训练集的scaler来转换验证集
        X2 = scaler.transform(X2)

        # Reshape data into the format accepted by the model
        X1 = X1.reshape(X1.shape[0], config.X_feature_len, 1)
        X1 = X1.astype('float64')
        X2 = X2.reshape(X2.shape[0], config.X_feature_len, 1)
        X2 = X2.astype('float64')

        # X1 = X1.reshape(X1.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X1 = X1.astype('float64')
        # X2 = X2.reshape(X2.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X2 = X2.astype('float64')

        # X2 = X2.reshape(X2.shape[0], config.lstm_seq_len, config.lstm_features_num)
        # X2 = X2.astype('float64')
        return X1, X2, Y, row_nums, matrix_files

