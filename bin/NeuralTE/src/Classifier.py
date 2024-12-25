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

from utils.show_util import showToolName, showTestParams
from keras.models import load_model
from configs import config
from configs import gpu_config
from DataProcessor import DataProcessor
from utils.evaluate_util import get_metrics, correct_using_minority
from utils.data_util import get_feature_len, get_gpu_config


class Classifier:
    def __init__(self, model_path):
        self.model = load_model(model_path)
        #self.model = load_model(model_path, compile=False)

    def predict(self, X, y, seq_names, data_path):
        # 预测概率
        y_pred = self.model.predict(X)
        accuracy, precision, recall, f1 = get_metrics(y_pred, y, seq_names, data_path)

        return accuracy, precision, recall, f1


def main():
    # showToolName()

    # 1.parse args
    describe_info = '########################## NeuralTE, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--data', required=True, metavar='data', help='Input fasta file used to predict, header format: seq_name\tlabel\tspecies_name, refer to "data/test.example.fa" for example.')
    parser.add_argument('--out_dir', required=True, metavar='output_dir', help='Output directory, store temporary files')
    parser.add_argument('--use_TSD', metavar='use_TSD', help='Whether to use TSD features, 1: true, 0: false. default = [ ' + str(config.use_TSD) + ' ]')
    parser.add_argument('--is_predict', metavar='is_predict', help='Enable prediction mode, 1: true, 0: false. default = [ ' + str(config.is_predict) + ' ]')

    parser.add_argument('--start_gpu_num', metavar='start_gpu_num', help='The starting index for using GPUs. default = [ ' + str(gpu_config.start_gpu_num) + ' ]')
    parser.add_argument('--use_gpu_num', metavar='use_gpu_num', help='Specifying the number of GPUs in use. default = [ ' + str(gpu_config.use_gpu_num) + ' ]')
    parser.add_argument('--keep_raw', metavar='keep_raw', help='Whether to retain the raw input sequence, 1: true, 0: false; only save species having TSDs. default = [ ' + str(config.keep_raw) + ' ]')
    parser.add_argument('--genome', metavar='genome', help='Genome path, use to search for TSDs')
    parser.add_argument('--species', metavar='species', help='Which species does the TE library to be classified come from.')
    parser.add_argument('--model_path', metavar='model_path', help='Input the path of trained model, absolute path.')
    parser.add_argument('--use_kmers', metavar='use_kmers', help='Whether to use kmers features, 1: true, 0: false. default = [ ' + str(config.use_kmers) + ' ]')
    parser.add_argument('--use_terminal', metavar='use_terminal', help='Whether to use LTR, TIR terminal features, 1: true, 0: false. default = [ ' + str(config.use_terminal) + ' ]')
    parser.add_argument('--use_minority', metavar='use_minority', help='Whether to use minority features, 1: true, 0: false. default = [ ' + str(config.use_minority) + ' ]')
    parser.add_argument('--use_domain', metavar='use_domain', help='Whether to use domain features, 1: true, 0: false. default = [ ' + str(config.use_domain) + ' ]')
    parser.add_argument('--use_ends', metavar='use_ends', help='Whether to use 5-bp terminal ends features, 1: true, 0: false. default = [ ' + str(config.use_ends) + ' ]')
    parser.add_argument('--is_wicker', metavar='is_wicker', help='Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [ ' + str(config.is_wicker) + ' ]')
    parser.add_argument('--is_plant', metavar='is_plant', help='Is the input genome of a plant? 0 represents non-plant, while 1 represents plant. default = [ ' + str(config.is_plant) + ' ]')
    parser.add_argument('--threads', metavar='thread_num', help='Input thread num, default = [ ' + str(config.threads) + ' ]')
    parser.add_argument('--internal_kmer_sizes', metavar='internal_kmer_sizes', help='The k-mer size used to convert internal sequences to k-mer frequency features, default = [ ' + str(config.internal_kmer_sizes) + ' MB ]')
    parser.add_argument('--terminal_kmer_sizes', metavar='terminal_kmer_sizes', help='The k-mer size used to convert terminal sequences to k-mer frequency features, default = [ ' + str(config.terminal_kmer_sizes) + ' ]')

    args = parser.parse_args()

    data_path = args.data
    out_dir = args.out_dir
    genome = args.genome
    species = args.species
    is_plant = args.is_plant
    model_path = args.model_path
    use_terminal = args.use_terminal
    use_kmers = args.use_kmers
    use_TSD = args.use_TSD
    use_minority = args.use_minority
    use_domain = args.use_domain
    use_ends = args.use_ends
    is_predict = args.is_predict

    start_gpu_num = args.start_gpu_num
    use_gpu_num = args.use_gpu_num
    keep_raw = args.keep_raw
    is_wicker = args.is_wicker
    threads = args.threads
    internal_kmer_sizes = args.internal_kmer_sizes
    terminal_kmer_sizes = args.terminal_kmer_sizes

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
    if is_predict is not None:
        config.is_predict = int(is_predict)

    if start_gpu_num is not None:
        gpu_config.start_gpu_num = int(start_gpu_num)
    if use_gpu_num is not None:
        gpu_config.use_gpu_num = int(use_gpu_num)
    if keep_raw is not None:
        config.keep_raw = int(keep_raw)
    if is_wicker is not None:
        config.is_wicker = int(is_wicker)
    if is_plant is not None:
        config.is_plant = int(is_plant)
    if threads is not None:
        config.threads = int(threads)
    if internal_kmer_sizes is not None:
        config.internal_kmer_sizes = json.loads(internal_kmer_sizes)
    if terminal_kmer_sizes is not None:
        config.terminal_kmer_sizes = json.loads(terminal_kmer_sizes)

    if genome is not None:
        os.makedirs(config.work_dir, exist_ok=True)
        genome_info_path = config.work_dir + '/genome.info'
        if str(genome).__contains__('genome.info'):
            os.system('cp ' + genome + ' ' + genome_info_path)
        else:
            if species is None:
                species = 'Unknown'
            with open(genome_info_path, 'w') as f_save:
                f_save.write('#Scientific Name\tGenome Path\tIs Plant\n')
                f_save.write(species+'\t'+genome+'\t'+str(config.is_plant)+'\n')

    params = {}
    params['data_path'] = data_path
    params['out_dir'] = out_dir
    params['model_path'] = model_path
    params['genome'] = genome
    params['species'] = species
    showTestParams(params)

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
    # print(X.shape, y.shape)
    endtime1 = time.time()
    dtime1 = endtime1 - starttime1
    print("Running time of DataProcessor: %.8s s" % (dtime1))

    starttime2 = time.time()
    # Instantiate the Classifier class
    classifier = Classifier(model_path=model_path)
    # start prediction
    accuracy, precision, recall, f1 = classifier.predict(X, y, seq_names, data_path)
    endtime2 = time.time()
    dtime2 = endtime2 - starttime2
    print("Running time of model Classifier: %.8s s" % (dtime2))

    if config.use_minority:
        starttime3 = time.time()
        correct_using_minority(data_path, threads)
        endtime3 = time.time()
        dtime3 = endtime3 - starttime3
        print("Running time of model Classifier: %.8s s" % (dtime3))

    endtime3 = time.time()
    dtime = endtime3 - starttime1
    print("Running time of total NeuralTE-Classifier: %.8s s" % (dtime))

if __name__ == '__main__':
    main()