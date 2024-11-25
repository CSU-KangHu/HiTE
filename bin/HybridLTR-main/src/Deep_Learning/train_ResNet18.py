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

from configs import config, gpu_config

from typing import List, Tuple
from matplotlib.pyplot import imshow

import matplotlib.pyplot as plt
from tensorflow.keras.preprocessing import image
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Add, Dense, Activation, ZeroPadding2D, BatchNormalization, Flatten, Conv2D, AveragePooling2D, MaxPooling2D, GlobalMaxPooling2D
from tensorflow.keras.initializers import glorot_uniform
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.utils import to_categorical

from keras.callbacks import EarlyStopping
from keras.layers import Dense, Conv2D,  MaxPool2D, Flatten, GlobalAveragePooling2D,  BatchNormalization, Layer, Add
from keras.models import Sequential
from keras.models import Model
import tensorflow as tf

from utils.show_util import showToolName, showTestParams
from keras.models import load_model
from configs import config
from configs import gpu_config
from DataProcessor import DataProcessor
from utils.evaluate_util import get_metrics, correct_using_minority, get_metrics_by_label, get_metrics_v1
from utils.data_util import get_feature_len, get_gpu_config, find_files_recursively
from ResNet18 import ResNet18
import tensorflow as tf
from PIL import Image

import os
import numpy as np
import keras
from keras import layers
from tensorflow import data as tf_data
import matplotlib.pyplot as plt

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
    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
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





get_gpu_config(gpu_config.start_gpu_num, gpu_config.use_gpu_num)

# set input image parameters
image_size = config.image_shape
channels = 3
num_classes = 2

input_dir = '/home/hukang/left_LTR_real_dataset/clean_data/both_ends_frames_clean_expand_no_Drosophila'
in_folder = input_dir + '_images'
# convert_matrix2jpg(input_dir, in_folder)


img_width = image_size[0]
img_height = image_size[1]
batch_size = 64

validation_split = 0.1
test_split = 0.1

# 加载数据集
full_dataset = tf.keras.preprocessing.image_dataset_from_directory(
    in_folder,
    image_size=(img_width, img_height),
    batch_size=batch_size,
    validation_split=validation_split + test_split,
    label_mode='categorical',
    seed=123,
    subset="training"
)

# 划分训练集和验证集
num_examples = len(full_dataset)
num_validation_examples = int(num_examples * validation_split)
num_test_examples = int(num_examples * test_split)

val_ds = full_dataset.take(num_validation_examples)
test_ds = full_dataset.skip(num_validation_examples).take(num_test_examples)
train_ds = full_dataset.skip(num_validation_examples + num_test_examples)

# train_ds = tf.keras.preprocessing.image_dataset_from_directory(
#     in_folder,
#     validation_split = 0.2,
#     subset="training",
#     label_mode='categorical', # default mode is 'int' label, but we want one-hot encoded labels (e.g. for categorical_crossentropy loss)
#     seed=123,
#     image_size=(img_width, img_height),
#     batch_size=batch_size
# )
#
# val_ds = tf.keras.preprocessing.image_dataset_from_directory(
#     in_folder,
#     validation_split=0.2,
#     subset="validation",
#     label_mode='categorical',
#     seed=123,
#     image_size=(img_width, img_height),
#     batch_size=batch_size
# )
print("Dataset shape:", train_ds)

# plt.figure(figsize=(10, 10))
# for images, labels in train_ds.take(1):
#     for i in range(9):
#         ax = plt.subplot(3, 3, i + 1)
#         plt.imshow(np.array(images[i]).astype("uint8"))
#         plt.show()
#         # plt.title(labels[i])
#         # plt.axis("off")

# class_names = train_ds.class_names
# print(class_names)

# use keras functionality for adding a rescaling layer
normalization_layer = tf.keras.layers.experimental.preprocessing.Rescaling(1./255)

# rescale training and validation sets
norm_train_ds = train_ds.map(lambda x, y: (normalization_layer(x), y))
norm_val_ds = val_ds.map(lambda x, y: (normalization_layer(x), y))
norm_test_ds = test_ds.map(lambda x, y: (normalization_layer(x), y))


AUTOTUNE = tf.data.AUTOTUNE

norm_train_ds = norm_train_ds.cache().prefetch(buffer_size=AUTOTUNE)
norm_val_ds = norm_val_ds.cache().prefetch(buffer_size=AUTOTUNE)
norm_test_ds = norm_test_ds.cache().prefetch(buffer_size=AUTOTUNE)

reduce_lr = tf.keras.callbacks.ReduceLROnPlateau(factor=0.1, patience=10, min_lr=0.00001)

callbacks = [
    # keras.callbacks.EarlyStopping(
    #     monitor="val_loss", # monitor validation loss (that is, the loss computed for the validation holdout)
    #     min_delta=1e-3, # "no longer improving" being defined as "an improvement lower than 1e-2"
    #     patience=10, # "no longer improving" being further defined as "for at least 10 consecutive epochs"
    #     verbose=1
    # ),
    reduce_lr
]

import time

start = time.time()

model = ResNet18(num_classes)
model.build(input_shape = (None, image_size[0], image_size[1], channels))
#use categorical_crossentropy since the label is one-hot encoded
from tensorflow.keras.optimizers import SGD
sgd = SGD(learning_rate=0.01, momentum=0.9, decay=1e-6, nesterov=True)
model.compile(optimizer = sgd,loss='categorical_crossentropy', metrics=["accuracy"])
model.summary()

history = model.fit(
    norm_train_ds,
    validation_data=norm_val_ds,
    epochs = 30,
    callbacks=callbacks
)

model_path = config.project_dir + '/models/ResNet18.h5'
model.save_weights(model_path)

print(history.history.keys())

plt.figure(figsize=(12, 10))

# summarize history for accuracy
plt.subplot(2, 1, 1)
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'val'], loc='upper left')

# summarize history for loss
plt.subplot(2, 1, 2)
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'val'], loc='upper left')
plt.show()

stop = time.time()

print(f'Training took: {(stop-start)/60} minutes')


model.load_weights(model_path)

y_pred = model.predict(norm_test_ds)
# 将预测结果转换为类别标签
y_pred = tf.argmax(y_pred, axis=1).numpy()

# 获取真实标签
y_test = []
for images, labels in norm_test_ds:
    y_test.extend(tf.argmax(labels, axis=1))
y_test = np.array(y_test)
accuracy, precision, recall, f1 = get_metrics_v1(y_pred, y_test)


