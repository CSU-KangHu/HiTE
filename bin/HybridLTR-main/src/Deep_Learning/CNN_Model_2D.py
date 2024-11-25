import os
import atexit
from keras.models import load_model, Model
from keras.layers import Input, Dense, Dropout, Flatten, Conv2D
from configs import config
from configs import gpu_config


import keras.backend as K
import tensorflow as tf


class CNN_Model_2D:
    def __init__(self, input_shape, output_class_num):
        self.input_shape = input_shape
        self.class_num = output_class_num

    def train(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def build_model(self, cnn_num_convs, cnn_filters_array):
        # Create a MirroredStrategy.
        strategy = tf.distribute.MirroredStrategy()
        # Open a strategy scope and create/restore the model
        with strategy.scope():
            # CNN model
            # input layer
            input_layer = Input(shape=self.input_shape)
            conv_input_layer = input_layer
            # Create multiple convolutional layers
            for i in range(cnn_num_convs):
                # Add convolutional layers
                conv = Conv2D(cnn_filters_array[i], config.cnn_kernel_sizes_array[i], activation='relu')(conv_input_layer)
                conv_input_layer = conv
            dropout1 = Dropout(0.5)(conv_input_layer)
            # Add flattening and fully connected layers
            flatten = Flatten()(dropout1)
            dense1 = Dense(128, activation='relu')(flatten)
            # dropout2 = Dropout(config.cnn_dropout)(dense1)
            # Output layer
            output_layer = Dense(int(config.class_num), activation='softmax')(dense1)
            # Build the model
            model = Model(inputs=input_layer, outputs=output_layer)
            # Compile the model
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            #model.compile(loss=[categorical_focal_loss(alpha=0.25, gamma=2)], optimizer='adam', metrics=['accuracy'])
            # Print model summary
            #model.summary()
        # atexit.register(strategy._extended._collective_ops._pool.close)  # type: ignore
        return model