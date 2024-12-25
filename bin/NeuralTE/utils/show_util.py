from configs import config
from configs import gpu_config

def showToolName():
    describe_image = '\n' + \
    '  _   _                      _ _______ ______ \n' + \
    ' | \ | |                    | |__   __|  ____| \n' + \
    ' |  \| | ___ _   _ _ __ __ _| |  | |  | |__ \n' + \
    ' | . ` |/ _ \ | | | \'__/ _` | |  | |  |  __| \n' + \
    ' | |\  |  __/ |_| | | | (_| | |  | |  | |____ \n' + \
    ' |_| \_|\___|\__,_|_|  \__,_|_|  |_|  |______| \n' + \
    'version ' + str(config.version_num) + '\n\n'
    print(describe_image)

def showTrainParams(params):
    print('\nParameters configuration\n'
        '====================================System settings========================================\n'
        '  [Setting] Input data used to train model = [ ' + str(params['data_path']) + ' ]\n'
        '  [Setting] Output directory = [ ' + str(params['out_dir']) + ' ]\n'
        '  [Setting] Genome path = [ ' + str(params['genome']) + ' ]\n'
        '  [Setting] Whether to use TSD features = [ ' + str(config.use_TSD) + ' ]\n'
        '  [Setting] Enable train mode = [ ' + str(config.is_train) + ' ]\n'
        '  [Setting] Enable prediction mode = [ ' + str(config.is_predict) + ' ]\n' 
        '  [Setting] Whether to retain the raw input sequence = [ ' + str(config.keep_raw) + ' ]\n'  
        '  [Setting] Whether to use kmers features = [ ' + str(config.use_kmers) + ' ]\n'                                                                                                                                               
        '  [Setting] Whether to use LTR, TIR terminal features = [ ' + str(config.use_terminal) + ' ]\n'
        '  [Setting] Whether to use minority features = [ ' + str(config.use_minority) + ' ]\n'
        '  [Setting] Whether to use domain features = [ ' + str(config.use_domain) + ' ]\n'
        '  [Setting] Whether to use 5-bp terminal ends features = [ ' + str(config.use_ends) + ' ]\n'
        '  [Setting] Input thread num = [ ' + str(config.threads) + ' ]\n'
        '  [Setting] The k-mer size used to convert internal sequences to k-mer frequency features = [ ' + str(config.internal_kmer_sizes) + ' ]\n'
        '  [Setting] The k-mer size used to convert terminal sequences to k-mer frequency features = [ ' + str(config.terminal_kmer_sizes) + ' ]\n'
        '  [Setting] The number of CNN convolutional layers = [ ' + str(config.cnn_num_convs) + ' ]\n'
        '  [Setting] The number of filters in each CNN convolutional layer = [ ' + str(config.cnn_filters_array) + ' ]\n'
        '  [Setting] The kernel size in each of CNN convolutional layer = [ ' + str(config.cnn_kernel_sizes_array) + ' ]\n'
        '  [Setting] The threshold of CNN Dropout = [ ' + str(config.cnn_dropout) + ' ]\n'
        '  [Setting] The batch size in training model = [ ' + str(config.batch_size) + ' ]\n'
        '  [Setting] The number of epochs in training model = [ ' + str(config.epochs) + ' ]\n'
        '  [Setting] Whether to use breakpoint training = [ ' + str(config.use_checkpoint) + ' ]'
        '  [Setting] The starting index for using GPUs = [ ' + str(gpu_config.start_gpu_num) + ' ]\n'
        '  [Setting] The number of GPUs in use = [ ' + str(gpu_config.use_gpu_num) + ' ]\n'
    )

def showTestParams(params):
    print('\nParameters configuration\n'
          '====================================System settings========================================\n'
          '  [Setting] Input data to be classified = [ ' + str(params['data_path']) + ' ]\n'
          '  [Setting] Output directory = [ ' + str(params['out_dir']) + ' ]\n'
          '  [Setting] Genome path = [ ' + str(params['genome']) + ' ]\n'   
          '  [Setting] species = [ ' + str(params['species']) + ' ]\n'                                                            
          '  [Setting] Input the path of trained model, absolute path = [ ' + str(params['model_path']) + ' ]\n'
          '  [Setting] Whether to use TSD features = [ ' + str(config.use_TSD) + ' ]\n'
          '  [Setting] Enable prediction mode = [ ' + str(config.is_predict) + ' ]\n'                                                                                                 
          '  [Setting] Whether to retain the raw input sequence = [ ' + str(config.keep_raw) + ' ]\n'  
          '  [Setting] Whether to use kmers features = [ ' + str(config.use_kmers) + ' ]\n'                                                                                                                                               
          '  [Setting] Whether to use LTR, TIR terminal features = [ ' + str(config.use_terminal) + ' ]\n'
          '  [Setting] Whether to use minority features = [ ' + str(config.use_minority) + ' ]\n'
          '  [Setting] Whether to use domain features = [ ' + str(config.use_domain) + ' ]\n'
          '  [Setting] Whether to use 5-bp terminal ends features = [ ' + str(config.use_ends) + ' ]\n' 
          '  [Setting] Use Wicker or RepeatMasker classification labels = [ ' + str(config.is_wicker) + ' ]\n' 
          '  [Setting] Is the input genome of a plant = [ ' + str(config.is_plant) + ' ]\n'                                                                                              
          '  [Setting] Input thread num = [ ' + str(config.threads) + ' ]'
          '  [Setting] The k-mer size used to convert internal sequences to k-mer frequency features = [ ' + str(config.internal_kmer_sizes) + ' ]\n'
          '  [Setting] The k-mer size used to convert terminal sequences to k-mer frequency features = [ ' + str(config.terminal_kmer_sizes) + ' ]\n'
          '  [Setting] The starting index for using GPUs = [ ' + str(gpu_config.start_gpu_num) + ' ]\n'
          '  [Setting] The number of GPUs in use = [ ' + str(gpu_config.use_gpu_num) + ' ]\n'
    )