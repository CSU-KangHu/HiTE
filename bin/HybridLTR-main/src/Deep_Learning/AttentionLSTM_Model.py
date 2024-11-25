import tensorflow as tf
from tensorflow.keras import layers, Model, initializers
from keras.models import load_model, Model
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense
from keras.layers import Input, Dense, Dropout, Flatten, Conv1D

import keras
import keras.backend as K
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Reshape, Flatten, LSTM, Dense, Dropout, Embedding, Bidirectional, GRU
from keras import initializers, regularizers, optimizers

from layers import AttentionWithContext, Addition

from configs import config


class AttentionLSTM:
    def __init__(self, seq_len, input_features_num, output_class_num):
        self.seq_len = seq_len
        self.num_features = input_features_num
        self.class_num = output_class_num

    def train(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def build_model(self, num_lstm_layers, hidden_size, is_attention, is_bidirectional):
        # Create a MirroredStrategy.
        strategy = tf.distribute.MirroredStrategy()
        # Open a strategy scope and create/restore the model
        with strategy.scope():
            num_classes = 2

            model = Sequential()
            input_layer = Input(shape=(self.seq_len, self.num_features))
            # input_layer = tf.squeeze(input_layer)
            model.add(input_layer)

            for i in range(num_lstm_layers):
                return_sequences = is_attention or (num_lstm_layers > 1 and i < num_lstm_layers - 1)
                print(return_sequences)
                if is_bidirectional:
                    model.add(Bidirectional(LSTM(hidden_size, return_sequences=return_sequences,
                                                 kernel_initializer=initializers.glorot_normal(seed=777),
                                                 bias_initializer='zeros')))
                else:
                    model.add(LSTM(hidden_size, return_sequences=return_sequences,
                                   kernel_initializer=initializers.glorot_normal(seed=777), bias_initializer='zeros'))

                if is_attention:
                    model.add(AttentionWithContext())
                    model.add(Addition())

            model.add(Dense(num_classes, activation='softmax', kernel_initializer=initializers.glorot_normal(seed=777),
                            bias_initializer='zeros'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            model.summary()

        return model