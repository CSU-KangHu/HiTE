import argparse
import json
import os
import sys
import time

from sklearn.preprocessing import StandardScaler
from tensorflow import keras



current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "../..")
sys.path.append(configs_folder)

import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
from keras.utils import np_utils
from configs import config, gpu_config
from Hybrid_Model import Hybrid_Model
from CNN_Model import CNN_Model
from DataProcessor import DataProcessor
from utils.evaluate_util import get_metrics, get_metrics_v1
from utils.show_util import showToolName, showTrainParams
from utils.data_util import get_gpu_config


class CrossValidator:
    def __init__(self, num_folds=5):
        self.num_folds = num_folds
        self.sss = StratifiedShuffleSplit(n_splits=num_folds, test_size=0.2, random_state=42)

    def evaluate(self, X, y, row_nums, matrix_files):
        accuracy_array = []
        precision_array = []
        recall_array = []
        f1_array = []
        copy_accuracy_array = []
        copy_precision_array = []
        copy_recall_array = []
        copy_f1_array = []
        # Loop through each K-fold
        for fold, (train_index, test_index) in enumerate(self.sss.split(X, y)):
            X_train, X_test = X[train_index], X[test_index]
            # X1_train, X1_test = X1[train_index], X1[test_index]
            # X2_train, X2_test = X2[train_index], X2[test_index]
            y_train, y_test = y[train_index], y[test_index]

            row_nums_train, row_nums_test = row_nums[train_index], row_nums[test_index]
            matrix_files_train, matrix_files_test= matrix_files[train_index], matrix_files[test_index]

            y_train_one_hot = np_utils.to_categorical(y_train, int(config.class_num))
            y_train_one_hot = np.array(y_train_one_hot)
            y_test_one_hot = np_utils.to_categorical(y_test, int(config.class_num))
            y_test_one_hot = np.array(y_test_one_hot)

            # 从训练集中再次划分出20%作为内部验证集
            split_idx = int(0.8 * len(X_train))
            X_train_inner, X_val_inner = X_train[:split_idx], X_train[split_idx:]
            y_train_inner, y_val_inner = y_train_one_hot[:split_idx], y_train_one_hot[split_idx:]


            # 可以在这里添加额外的逻辑，比如评估模型在外部验证集上的性能
            # val_loss, val_acc = model.evaluate(X_val, y_val)
            # print(f'Validation loss: {val_loss}, Validation accuracy: {val_acc}')


            # hybrid_model = Hybrid_Model(config.class_num)
            # model = hybrid_model.build_model()

            cnn_model = CNN_Model(config.X_feature_len, config.class_num)
            model = cnn_model.build_model(config.cnn_num_convs, config.cnn_filters_array)

            # cnn_model = CNN_Model_2D(config.input_shape, config.class_num)
            # model = cnn_model.build_model(config.cnn_num_convs, config.cnn_filters_array)


            # # CNN 模型
            # cnn_model = CNN_Model(config.image_shape, config.class_num)
            # model = cnn_model.build_model(config.cnn_num_convs, config.cnn_filters_array)

            # # LSTM 模型
            # lstm_model = LSTM_Model(config.lstm_seq_len, config.lstm_features_num, config.class_num)
            # model = lstm_model.build_model(config.lstm_num_layers, config.lstm_hidden_size, config.lstm_dropout)

            # # AttentionLSTM 模型
            # lstm_model = AttentionLSTM(config.lstm_seq_len, config.lstm_features_num, config.class_num)
            # model = lstm_model.build_model(config.lstm_num_layers, config.lstm_hidden_size, True, True)

            reduce_lr = keras.callbacks.ReduceLROnPlateau(factor=0.1, patience=3, verbose=1)

            # Train the model
            model.fit(X_train_inner, y_train_inner, batch_size=config.batch_size,
                      epochs=config.epochs, validation_data=(X_val_inner, y_val_inner),
                      verbose=1, callbacks=[reduce_lr])
            # Save the model
            model_path = config.project_dir + '/models/' + f'model_fold_{fold}.h5'
            model.save(model_path)

            # Predict probabilities
            y_pred = model.predict(X_test)
            accuracy, precision, recall, f1, copy_accuracy, copy_precision, copy_recall, copy_f1 = get_metrics(y_pred, y_test, row_nums_test, matrix_files_test)
            print("Fold:", fold)
            accuracy_array.append(accuracy)
            precision_array.append(precision)
            recall_array.append(recall)
            f1_array.append(f1)
            copy_accuracy_array.append(copy_accuracy)
            copy_precision_array.append(copy_precision)
            copy_recall_array.append(copy_recall)
            copy_f1_array.append(copy_f1)
        accuracies = np.array(accuracy_array)
        precisions = np.array(precision_array)
        recalls = np.array(recall_array)
        f1s = np.array(f1_array)
        copy_accuracies = np.array(copy_accuracy_array)
        copy_precisions = np.array(copy_precision_array)
        copy_recalls = np.array(copy_recall_array)
        copy_f1s = np.array(copy_f1_array)

        # Calculate mean and sample standard deviation
        accuracy_mean = round(np.mean(accuracies), 4)
        accuracy_stdv = round(np.std(accuracies, ddof=1), 4)  # Use ddof=1 to calculate sample standard deviation
        precision_mean = round(np.mean(precisions), 4)
        precision_stdv = round(np.std(precisions, ddof=1), 4)
        recall_mean = round(np.mean(recalls), 4)
        recall_stdv = round(np.std(recalls, ddof=1), 4)
        f1_mean = round(np.mean(f1s), 4)
        f1_stdv = round(np.std(f1s, ddof=1), 4)

        # Calculate mean and sample standard deviation
        copy_accuracy_mean = round(np.mean(copy_accuracies), 4)
        copy_accuracy_stdv = round(np.std(copy_accuracies, ddof=1), 4)  # Use ddof=1 to calculate sample standard deviation
        copy_precision_mean = round(np.mean(copy_precisions), 4)
        copy_precision_stdv = round(np.std(copy_precisions, ddof=1), 4)
        copy_recall_mean = round(np.mean(copy_recalls), 4)
        copy_recall_stdv = round(np.std(copy_recalls, ddof=1), 4)
        copy_f1_mean = round(np.mean(copy_f1s), 4)
        copy_f1_stdv = round(np.std(copy_f1s, ddof=1), 4)
        return [accuracy_mean, precision_mean, recall_mean, f1_mean], [accuracy_stdv, precision_stdv, recall_stdv, f1_stdv], [copy_accuracy_mean, copy_precision_mean, copy_recall_mean, copy_f1_mean], [copy_accuracy_stdv, copy_precision_stdv, copy_recall_stdv, copy_f1_stdv]


    def evaluate_hybrid(self, X1, X2, y, row_nums, matrix_files):
        accuracy_array = []
        precision_array = []
        recall_array = []
        f1_array = []
        copy_accuracy_array = []
        copy_precision_array = []
        copy_recall_array = []
        copy_f1_array = []

        # Loop through each K-fold
        for fold, (train_index, test_index) in enumerate(self.sss.split(X1, y)):
            X1_train, X1_test = X1[train_index], X1[test_index]
            X2_train, X2_test = X2[train_index], X2[test_index]
            y_train, y_test = y[train_index], y[test_index]

            row_nums_train, row_nums_test = row_nums[train_index], row_nums[test_index]
            matrix_files_train, matrix_files_test= matrix_files[train_index], matrix_files[test_index]

            y_train_one_hot = np_utils.to_categorical(y_train, int(config.class_num))
            y_train_one_hot = np.array(y_train_one_hot)
            y_test_one_hot = np_utils.to_categorical(y_test, int(config.class_num))
            y_test_one_hot = np.array(y_test_one_hot)

            # 从训练集中再次划分出20%作为内部验证集
            split_idx = int(0.8 * len(X1_train))
            X1_train_inner, X1_val_inner = X1_train[:split_idx], X1_train[split_idx:]
            X2_train_inner, X2_val_inner = X2_train[:split_idx], X2_train[split_idx:]
            y_train_inner, y_val_inner = y_train_one_hot[:split_idx], y_train_one_hot[split_idx:]


            hybrid_model = Hybrid_Model(config.class_num)
            model = hybrid_model.build_model()

            reduce_lr = keras.callbacks.ReduceLROnPlateau(factor=0.1, patience=3, verbose=1)

            # Train the model
            model.fit([X1_train_inner, X2_train_inner], y_train_inner, batch_size=config.batch_size,
                      epochs=config.epochs, validation_data=([X1_val_inner, X2_val_inner], y_val_inner),
                      verbose=1, callbacks=[reduce_lr])
            # Save the model
            model_path = config.project_dir + '/models/' + f'model_fold_{fold}.h5'
            model.save(model_path)

            # Predict probabilities
            y_pred = model.predict([X1_test, X2_test])
            accuracy, precision, recall, f1, copy_accuracy, copy_precision, copy_recall, copy_f1 = get_metrics(y_pred, y_test, row_nums_test, matrix_files_test)
            print("Fold:", fold)
            accuracy_array.append(accuracy)
            precision_array.append(precision)
            recall_array.append(recall)
            f1_array.append(f1)
            copy_accuracy_array.append(copy_accuracy)
            copy_precision_array.append(copy_precision)
            copy_recall_array.append(copy_recall)
            copy_f1_array.append(copy_f1)
        accuracies = np.array(accuracy_array)
        precisions = np.array(precision_array)
        recalls = np.array(recall_array)
        f1s = np.array(f1_array)
        copy_accuracies = np.array(copy_accuracy_array)
        copy_precisions = np.array(copy_precision_array)
        copy_recalls = np.array(copy_recall_array)
        copy_f1s = np.array(copy_f1_array)

        # Calculate mean and sample standard deviation
        accuracy_mean = round(np.mean(accuracies), 4)
        accuracy_stdv = round(np.std(accuracies, ddof=1), 4)  # Use ddof=1 to calculate sample standard deviation
        precision_mean = round(np.mean(precisions), 4)
        precision_stdv = round(np.std(precisions, ddof=1), 4)
        recall_mean = round(np.mean(recalls), 4)
        recall_stdv = round(np.std(recalls, ddof=1), 4)
        f1_mean = round(np.mean(f1s), 4)
        f1_stdv = round(np.std(f1s, ddof=1), 4)

        # Calculate mean and sample standard deviation
        copy_accuracy_mean = round(np.mean(copy_accuracies), 4)
        copy_accuracy_stdv = round(np.std(copy_accuracies, ddof=1), 4)  # Use ddof=1 to calculate sample standard deviation
        copy_precision_mean = round(np.mean(copy_precisions), 4)
        copy_precision_stdv = round(np.std(copy_precisions, ddof=1), 4)
        copy_recall_mean = round(np.mean(copy_recalls), 4)
        copy_recall_stdv = round(np.std(copy_recalls, ddof=1), 4)
        copy_f1_mean = round(np.mean(copy_f1s), 4)
        copy_f1_stdv = round(np.std(copy_f1s, ddof=1), 4)
        return [accuracy_mean, precision_mean, recall_mean, f1_mean], [accuracy_stdv, precision_stdv, recall_stdv, f1_stdv], [copy_accuracy_mean, copy_precision_mean, copy_recall_mean, copy_f1_mean], [copy_accuracy_stdv, copy_precision_stdv, copy_recall_stdv, copy_f1_stdv]

    def evaluate_hybrid_v1(self, X1, X2, y):
        accuracy_array = []
        precision_array = []
        recall_array = []
        f1_array = []
        # Loop through each K-fold
        for fold, (train_index, test_index) in enumerate(self.sss.split(X1, y)):
            X1_train, X1_test = X1[train_index], X1[test_index]
            X2_train, X2_test = X2[train_index], X2[test_index]
            y_train, y_test = y[train_index], y[test_index]

            y_train_one_hot = np_utils.to_categorical(y_train, int(config.class_num))
            y_train_one_hot = np.array(y_train_one_hot)
            y_test_one_hot = np_utils.to_categorical(y_test, int(config.class_num))
            y_test_one_hot = np.array(y_test_one_hot)

            # 从训练集中再次划分出20%作为内部验证集
            split_idx = int(0.8 * len(X1_train))
            X1_train_inner, X1_val_inner = X1_train[:split_idx], X1_train[split_idx:]
            X2_train_inner, X2_val_inner = X2_train[:split_idx], X2_train[split_idx:]
            y_train_inner, y_val_inner = y_train_one_hot[:split_idx], y_train_one_hot[split_idx:]

            hybrid_model = Hybrid_Model(config.class_num)
            model = hybrid_model.build_model()

            reduce_lr = keras.callbacks.ReduceLROnPlateau(factor=0.1, patience=3, verbose=1)

            # Train the model
            model.fit([X1_train_inner, X2_train_inner], y_train_inner, batch_size=config.batch_size,
                      epochs=config.epochs, validation_data=([X1_val_inner, X2_val_inner], y_val_inner),
                      verbose=1, callbacks=[reduce_lr])
            # Save the model
            model_path = config.project_dir + '/models/' + f'model_fold_{fold}.h5'
            model.save(model_path)

            # Predict probabilities
            y_pred = model.predict([X1_test, X2_test])
            accuracy, precision, recall, f1 = get_metrics_v1(y_pred, y_test)
            print("Fold:", fold)
            accuracy_array.append(accuracy)
            precision_array.append(precision)
            recall_array.append(recall)
            f1_array.append(f1)
        accuracies = np.array(accuracy_array)
        precisions = np.array(precision_array)
        recalls = np.array(recall_array)
        f1s = np.array(f1_array)

        # Calculate mean and sample standard deviation
        accuracy_mean = round(np.mean(accuracies), 4)
        accuracy_stdv = round(np.std(accuracies, ddof=1), 4)  # Use ddof=1 to calculate sample standard deviation
        precision_mean = round(np.mean(precisions), 4)
        precision_stdv = round(np.std(precisions, ddof=1), 4)
        recall_mean = round(np.mean(recalls), 4)
        recall_stdv = round(np.std(recalls, ddof=1), 4)
        f1_mean = round(np.mean(f1s), 4)
        f1_stdv = round(np.std(f1s, ddof=1), 4)

        return [accuracy_mean, precision_mean, recall_mean, f1_mean], [accuracy_stdv, precision_stdv, recall_stdv, f1_stdv]



def main():
    showToolName()

    # 1.parse args
    describe_info = '########################## NeuralTE-CrossValidator, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data_dir', required=True, metavar='data_dir',
                        help='Input directory contains positive and negative data')
    parser.add_argument('--out_dir', required=True, metavar='output_dir',
                        help='Output directory, store temporary files')

    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')

    parser.add_argument('--threads', metavar='thread_num', help='Input thread num, default = [ ' + str(config.threads) + ' ]')

    parser.add_argument('--cnn_num_convs', metavar='cnn_num_convs', help='The number of CNN convolutional layers. default = [ ' + str(config.cnn_num_convs) + ' ]')
    parser.add_argument('--cnn_filters_array', metavar='cnn_filters_array', help='The number of filters in each CNN convolutional layer. default = [ ' + str(config.cnn_filters_array) + ' ]')
    parser.add_argument('--cnn_kernel_sizes_array', metavar='cnn_kernel_sizes_array', help='The kernel size in each of CNN convolutional layer. default = [ ' + str(config.cnn_kernel_sizes_array) + ' ]')
    parser.add_argument('--cnn_dropout', metavar='cnn_dropout', help='The threshold of CNN Dropout. default = [ ' + str(config.cnn_dropout) + ' ]')
    parser.add_argument('--batch_size', metavar='batch_size', help='The batch size in training model. default = [ ' + str(config.batch_size) + ' ]')
    parser.add_argument('--epochs', metavar='epochs', help='The number of epochs in training model. default = [ ' + str(config.epochs) + ' ]')
    parser.add_argument('--use_checkpoint', metavar='use_checkpoint', help='Whether to use breakpoint training. 1: true, 0: false. The model will continue training from the last failed parameters to avoid training from head. default = [ ' + str(config.use_checkpoint) + ' ]')

    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    threads = args.threads
    start_gpu_num = args.start_gpu_num
    use_gpu_num = args.use_gpu_num

    cnn_num_convs = args.cnn_num_convs
    cnn_filters_array = args.cnn_filters_array
    cnn_kernel_sizes_array = args.cnn_kernel_sizes_array
    cnn_dropout = args.cnn_dropout
    batch_size = args.batch_size
    epochs = args.epochs
    use_checkpoint = args.use_checkpoint

    if out_dir is not None:
        config.work_dir = out_dir

    if threads is not None:
        config.threads = int(threads)
    if start_gpu_num is not None:
        gpu_config.start_gpu_num = int(start_gpu_num)
    if use_gpu_num is not None:
        gpu_config.use_gpu_num = int(use_gpu_num)

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
    X_feature_len = 0
    for kmer_size in config.kmer_sizes:
        X_feature_len += pow(5, kmer_size) - 1
    # X_feature_len += 5 * 100
    config.X_feature_len = X_feature_len

    # X_feature_len = 0
    # for kmer_size in config.kmer_sizes:
    #     X_feature_len += (100 - kmer_size + 1) * 1
    # config.X_feature_len = X_feature_len

    # lstm_seq_len = 0
    # for kmer_size in config.kmer_sizes:
    #     lstm_seq_len += 100 - kmer_size + 1
    # config.lstm_seq_len = lstm_seq_len
    # config.lstm_features_num = 2

    # reload GPU config
    get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

    # Instantiate the DataProcessor class
    data_processor = DataProcessor()
    # Load the data
    # X, y, row_nums, matrix_files = data_processor.load_data(positive_dir, negative_dir)
    # print(X.shape, y.shape)

    # Load the data
    X1, X2, y, row_nums, matrix_files = data_processor.load_data_hybrid(positive_dir, negative_dir)
    print(X1.shape, X2.shape, y.shape)



    # # 保存数组到文件
    # np.save(out_dir + '/X1.npy', X1)
    # np.save(out_dir + '/X2.npy', X2)
    # np.save(out_dir + '/y.npy', y)

    # X1 = np.load(out_dir + '/X1.npy')
    # X2 = np.load(out_dir + '/X2.npy')
    # y = np.load(out_dir + '/y.npy')
    #
    # X1 = X1[:, :, :1]
    # X2 = X2[:, :, :1]
    # print(X1.shape, X2.shape, y.shape)
    #
    # X1 = X1.reshape(X1.shape[0], config.X_feature_len, 1)
    # X1 = X1.astype('float64')
    # X2 = X2.reshape(X2.shape[0], config.X_feature_len, 1)
    # X2 = X2.astype('float64')

    starttime2 = time.time()

    # Instantiate the CrossValidator class
    validator = CrossValidator(num_folds=5)
    # Perform cross-validation
    # means, stdvs, copy_means, copy_stdvs = validator.evaluate(X, y, row_nums, matrix_files)
    means, stdvs, copy_means, copy_stdvs = validator.evaluate_hybrid(X1, X2, y, row_nums, matrix_files)
    # means, stdvs = validator.evaluate_hybrid_v1(X1, X2, y)
    print('accuracy, precision, recall, f1:')
    print("Mean array:", means)
    print("stdv array:", stdvs)

    print('copy accuracy, precision, recall, f1:')
    print("Mean array:", copy_means)
    print("stdv array:", copy_stdvs)

    endtime2 = time.time()
    dtime2 = endtime2 - starttime2
    print("Running time of model Trainer: %.8s s" % (dtime2))


if __name__ == '__main__':
    main()