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
        '  [Setting] Input data used to train model = [ ' + str(params['data_dir']) + ' ]\n'
        '  [Setting] Output directory = [ ' + str(params['out_dir']) + ' ]\n'
        '  [Setting] Input thread num = [ ' + str(config.threads) + ' ]\n'
        '  [Setting] The number of CNN convolutional layers = [ ' + str(config.cnn_num_convs) + ' ]\n'
        '  [Setting] The number of filters in each CNN convolutional layer = [ ' + str(config.cnn_filters_array) + ' ]\n'
        '  [Setting] The kernel size in each of CNN convolutional layer = [ ' + str(config.cnn_kernel_sizes_array) + ' ]\n'
        '  [Setting] The threshold of CNN Dropout = [ ' + str(config.cnn_dropout) + ' ]\n'
        '  [Setting] The batch size in training model = [ ' + str(config.batch_size) + ' ]\n'
        '  [Setting] The number of epochs in training model = [ ' + str(config.epochs) + ' ]\n'
        '  [Setting] Whether to use breakpoint training = [ ' + str(config.use_checkpoint) + ' ]\n'
        '  [Setting] The starting index for using GPUs = [ ' + str(gpu_config.start_gpu_num) + ' ]\n'
        '  [Setting] The number of GPUs in use = [ ' + str(gpu_config.use_gpu_num) + ' ]\n'
    )

def showTestParams(params):
    print('\nParameters configuration\n'
          '====================================System settings========================================\n'
          '  [Setting] Input data to be classified = [ ' + str(params['data_dir']) + ' ]\n'
          '  [Setting] Output directory = [ ' + str(params['out_dir']) + ' ]\n'                                                     
          '  [Setting] Input the path of trained model, absolute path = [ ' + str(params['model_path']) + ' ]\n'  
          '  [Setting] Input thread num = [ ' + str(config.threads) + ' ]\n'
          '  [Setting] The starting index for using GPUs = [ ' + str(gpu_config.start_gpu_num) + ' ]\n'
          '  [Setting] The number of GPUs in use = [ ' + str(gpu_config.use_gpu_num) + ' ]\n'
    )