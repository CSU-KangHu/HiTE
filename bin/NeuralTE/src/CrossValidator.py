import argparse
import os
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
from keras.utils import np_utils
from configs import config, gpu_config
from CNN_Model import CNN_Model
from DataProcessor import DataProcessor
from utils.evaluate_util import get_metrics
from utils.show_util import showToolName
from utils.data_util import get_feature_len, get_gpu_config


class CrossValidator:
    def __init__(self, num_folds=5):
        self.num_folds = num_folds
        self.sss = StratifiedShuffleSplit(n_splits=num_folds, test_size=0.2, random_state=42)

    def evaluate(self, X, y):
        accuracy_array = []
        precision_array = []
        recall_array = []
        f1_array = []
        # Loop through each K-fold
        for fold, (train_index, test_index) in enumerate(self.sss.split(X, y)):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            y_train_one_hot = np_utils.to_categorical(y_train, int(config.class_num))
            y_train_one_hot = np.array(y_train_one_hot)
            y_test_one_hot = np_utils.to_categorical(y_test, int(config.class_num))
            y_test_one_hot = np.array(y_test_one_hot)

            cnn_model = CNN_Model(config.X_feature_len, config.class_num)
            model = cnn_model.build_model(config.cnn_num_convs, config.cnn_filters_array)

            # Train the model
            model.fit(X_train, y_train_one_hot, batch_size=config.batch_size, epochs=config.epochs, verbose=1)
            # Save the model
            model_path = config.project_dir + '/models/' + f'model_fold_{fold}.h5'
            model.save(model_path)

            # Predict probabilities
            y_pred = model.predict(X_test)
            accuracy, precision, recall, f1 = get_metrics(y_pred, y_test, None, None)
            print("Fold:", fold)
            print("Accuracy:", accuracy)
            print("Precision:", precision)
            print("Recall:", recall)
            print("F1:", f1)
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
    parser.add_argument('--data', metavar='data', help='Input fasta file containing TSDs information used to CrossValidator model')
    parser.add_argument('--out_dir', metavar='output_dir', help='Output directory, store temporary files')
    parser.add_argument('--use_kmers', metavar='use_kmers', help='Whether to use kmers features, 1: true, 0: false. default = [ ' + str(config.use_kmers) + ' ]')
    parser.add_argument('--use_terminal', metavar='use_terminal', help='Whether to use LTR, TIR terminal features, 1: true, 0: false. default = [ ' + str(config.use_terminal) + ' ]')
    parser.add_argument('--use_TSD', metavar='use_TSD', help='Whether to use TSD features, 1: true, 0: false. default = [ ' + str(config.use_TSD) + ' ]')
    parser.add_argument('--use_domain', metavar='use_domain', help='Whether to use domain features, 1: true, 0: false. default = [ ' + str(config.use_domain) + ' ]')
    parser.add_argument('--use_ends', metavar='use_ends', help='Whether to use 5-bp terminal ends features, 1: true, 0: false. default = [ ' + str(config.use_ends) + ' ]')
    parser.add_argument('--is_predict', metavar='is_predict', help='Enable prediction mode, 1: true, 0: false. default = [ ' + str(config.is_predict) + ' ]')
    parser.add_argument('--threads', metavar='thread_num', help='Input thread num, default = [ ' + str(config.threads) + ' ]')
    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')

    args = parser.parse_args()

    data_path = args.data
    out_dir = args.out_dir
    use_kmers = args.use_kmers
    use_terminal = args.use_terminal
    use_TSD = args.use_TSD
    use_domain = args.use_domain
    use_ends = args.use_ends
    is_predict = args.is_predict
    threads = args.threads
    start_gpu_num = args.start_gpu_num
    use_gpu_num = args.use_gpu_num

    if out_dir is not None:
        config.work_dir = out_dir
    if use_terminal is not None:
        config.use_terminal = int(use_terminal)
    if use_kmers is not None:
        config.use_kmers = int(use_kmers)
    if use_TSD is not None:
        config.use_TSD = int(use_TSD)
    if use_domain is not None:
        config.use_domain = int(use_domain)
    if use_ends is not None:
        config.use_ends = int(use_ends)
    if is_predict is not None:
        config.is_predict = int(is_predict)
    if threads is not None:
        config.threads = int(threads)
    if start_gpu_num is not None:
        gpu_config.start_gpu_num = int(start_gpu_num)
    if use_gpu_num is not None:
        gpu_config.use_gpu_num = int(use_gpu_num)

    X_feature_len = get_feature_len()
    config.X_feature_len = X_feature_len

    # reload GPU config
    get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

    # Instantiate the DataProcessor class
    data_processor = DataProcessor()
    # Load the data
    # Make sure the header format for the following data is in Repbase format, i.e., 'TE_name  Superfamily Species', separated by '\t'
    cv_train_data_path = data_path  # Path to cross-validation training data
    X, y, seq_names, data_path = data_processor.load_data(config.internal_kmer_sizes, config.terminal_kmer_sizes, cv_train_data_path)
    print(X.shape, y.shape)

    # Instantiate the CrossValidator class
    validator = CrossValidator(num_folds=5)
    # Perform cross-validation
    means, stdvs = validator.evaluate(X, y)
    print('accuracy, precision, recall, f1:')
    print("Mean array:", means)
    print("stdv array:", stdvs)


if __name__ == '__main__':
    main()