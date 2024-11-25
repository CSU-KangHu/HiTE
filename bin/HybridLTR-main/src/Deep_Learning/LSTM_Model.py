import os
import atexit
from keras.models import load_model, Model
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, Bidirectional
from keras.layers import Input, Dense, Dropout, Flatten, Conv1D
from configs import config
from configs import gpu_config


import keras.backend as K
import tensorflow as tf


class LSTM_Model:
    def __init__(self, seq_len, input_features_num, output_class_num):
        self.seq_len = seq_len
        self.num_features = input_features_num
        self.class_num = output_class_num

    def train(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def build_model(self, num_lstm_layers, hidden_size, dropout_rate):
        # Create a MirroredStrategy.
        strategy = tf.distribute.MirroredStrategy()
        # Open a strategy scope and create/restore the model
        with strategy.scope():
            # LSTM model
            input_layer = Input(shape=(self.seq_len, self.num_features))
            lstm = input_layer
            for i in range(num_lstm_layers):
                if i == num_lstm_layers - 1:
                    lstm = Bidirectional(LSTM(units=hidden_size, dropout=dropout_rate, return_sequences=False))(lstm)
                else:
                    lstm = Bidirectional(LSTM(units=hidden_size, dropout=dropout_rate, return_sequences=True))(lstm)
            output_layer = Dense(int(config.class_num), activation='softmax')(lstm)
            model = Model(input_layer, output_layer)
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        # atexit.register(strategy._extended._collective_ops._pool.close)  # type: ignore
        return model