import os
import atexit
from keras.models import load_model, Model
from keras.layers import Input, Dense, Dropout, Flatten, Conv1D
from configs import config
from configs import gpu_config


import keras.backend as K
import tensorflow as tf


def categorical_focal_loss(gamma=2.0, alpha=0.25):
    """
    Implementation of Focal Loss from the paper in multiclass classification
    Formula:
        loss = -alpha*((1-p)^gamma)*log(p)
    Parameters:
        alpha -- the same as wighting factor in balanced cross entropy
        gamma -- focusing parameter for modulating factor (1-p)
    Default value:
        gamma -- 2.0 as mentioned in the paper
        alpha -- 0.25 as mentioned in the paper
    """

    def focal_loss(y_true, y_pred):
        # Define epsilon so that the backpropagation will not result in NaN
        # for 0 divisor case
        epsilon = K.epsilon()
        # Add the epsilon to prediction value
        # y_pred = y_pred + epsilon
        # Clip the prediction value
        y_pred = K.clip(y_pred, epsilon, 1.0 - epsilon)
        # Calculate cross entropy
        cross_entropy = -y_true * K.log(y_pred)
        # Calculate weight that consists of  modulating factor and weighting factor
        weight = alpha * y_true * K.pow((1 - y_pred), gamma)
        # Calculate focal loss
        loss = weight * cross_entropy
        # Sum the losses in mini_batch
        loss = K.sum(loss, axis=1)
        return loss

    return focal_loss


def binary_focal_loss(gamma=2.0, alpha=0.25):
    """
    Implementation of Focal Loss from the paper in multiclass classification
    Formula:
        loss = -alpha_t*((1-p_t)^gamma)*log(p_t)

        p_t = y_pred, if y_true = 1
        p_t = 1-y_pred, otherwise

        alpha_t = alpha, if y_true=1
        alpha_t = 1-alpha, otherwise

        cross_entropy = -log(p_t)
    Parameters:
        alpha -- the same as wighting factor in balanced cross entropy
        gamma -- focusing parameter for modulating factor (1-p)
    Default value:
        gamma -- 2.0 as mentioned in the paper
        alpha -- 0.25 as mentioned in the paper
    """

    def focal_loss(y_true, y_pred):
        # Define epsilon so that the backpropagation will not result in NaN
        # for 0 divisor case
        epsilon = K.epsilon()
        # Add the epsilon to prediction value
        # y_pred = y_pred + epsilon
        # Clip the prediciton value
        y_pred = K.clip(y_pred, epsilon, 1.0 - epsilon)
        # Calculate p_t
        p_t = tf.where(K.equal(y_true, 1), y_pred, 1 - y_pred)
        # Calculate alpha_t
        alpha_factor = K.ones_like(y_true) * alpha
        alpha_t = tf.where(K.equal(y_true, 1), alpha_factor, 1 - alpha_factor)
        # Calculate cross entropy
        cross_entropy = -K.log(p_t)
        weight = alpha_t * K.pow((1 - p_t), gamma)
        # Calculate focal loss
        loss = weight * cross_entropy
        # Sum the losses in mini_batch
        loss = K.sum(loss, axis=1)
        return loss

    return focal_loss

class CNN_Model:
    def __init__(self, input_features_num, output_class_num):
        self.num_features = input_features_num
        self.class_num = output_class_num

    def train(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def build_model(self, cnn_num_convs, cnn_filters_array):
        # Prepare a directory to store all the checkpoints.
        checkpoint_dir = config.work_dir + "/ckpt"
        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)

        # construct model
        if config.use_checkpoint == 0:
            os.system('cd ' + gpu_config.checkpoint_dir + ' && rm -rf ckpt*')
        # Create a MirroredStrategy.
        strategy = tf.distribute.MirroredStrategy()
        # Open a strategy scope and create/restore the model
        with strategy.scope():
            # Either restore the latest model, or create a fresh one
            # if there is no checkpoint available.
            checkpoints = [gpu_config.checkpoint_dir + "/" + name for name in os.listdir(gpu_config.checkpoint_dir)]
            if checkpoints:
                latest_checkpoint = max(checkpoints, key=os.path.getctime)
                print("Restoring from", latest_checkpoint)
                return load_model(latest_checkpoint)
            print("Creating a new model")

            # CNN model
            # input layer
            input_layer = Input(shape=(self.num_features, 1))
            conv_input_layer = input_layer
            # Create multiple convolutional layers
            for i in range(cnn_num_convs):
                # Add convolutional layers
                conv = Conv1D(cnn_filters_array[i], config.cnn_kernel_sizes_array[i], activation='relu')(conv_input_layer)
                conv_input_layer = conv
            dropout1 = Dropout(0.5)(conv_input_layer)
            # Add flattening and fully connected layers
            flatten = Flatten()(dropout1)
            dense1 = Dense(128, activation='relu')(flatten)
            dropout2 = Dropout(config.cnn_dropout)(dense1)
            # Output layer
            output_layer = Dense(int(config.class_num), activation='softmax')(dropout2)
            # Build the model
            model = Model(inputs=input_layer, outputs=output_layer)
            # Compile the model
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            #model.compile(loss=[categorical_focal_loss(alpha=0.25, gamma=2)], optimizer='adam', metrics=['accuracy'])
            # Print model summary
            #model.summary()
        atexit.register(strategy._extended._collective_ops._pool.close)  # type: ignore
        return model