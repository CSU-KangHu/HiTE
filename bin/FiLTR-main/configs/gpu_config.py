# 1. Configuring multi-GPU usage
start_gpu_num = 0 # Start GPU numbering
use_gpu_num = 1 # Number of GPUs to use, start_gpu_num=0, use_gpu_num=2 indicates using GPU0, GPU1, totaling two GPUs
all_devices = ["/gpu:0", "/gpu:1", "/gpu:2", "/gpu:3", "/gpu:4", "/gpu:5", "/gpu:6", "/gpu:7"] # All GPU numbers on the machine. If your machine has more than 8 GPUs, add "/gpu:8", "/gpu:9" ... afterward.


# 2. The following parameters do not need modification
import os
import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR')

# # Set the GPUs you want to use
# gpus = tf.config.experimental.list_physical_devices('GPU')
# # For GPU memory growth
# for device in gpus:
#     tf.config.experimental.set_memory_growth(device, True)
# use_devices = all_devices[start_gpu_num: start_gpu_num + use_gpu_num]
# tf.config.experimental.set_visible_devices(gpus[start_gpu_num: start_gpu_num + use_gpu_num], 'GPU')
#
# # # Use CPU to train model
# # num_CPU_cores = 40
# # tf.config.threading.set_inter_op_parallelism_threads(num_CPU_cores)
# # tf.config.threading.set_intra_op_parallelism_threads(num_CPU_cores)
# # tf.config.experimental.set_visible_devices([], 'GPU')
#
# # Create a MirroredStrategy to use multiple GPUs
# strategy = tf.distribute.MirroredStrategy(devices=use_devices)
# # print('Number of devices: {}'.format(strategy.num_replicas_in_sync))


# # Prepare a directory to store all the checkpoints.
checkpoint_dir = "./ckpt"
# if not os.path.exists(checkpoint_dir):
#     os.makedirs(checkpoint_dir)