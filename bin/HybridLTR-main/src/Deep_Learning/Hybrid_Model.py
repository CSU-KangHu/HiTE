import os

from keras.models import load_model, Model
from keras.layers import Input, Dense, Dropout, Flatten, Conv1D, Conv2D, LSTM,\
    Embedding, GlobalAveragePooling1D, concatenate, \
    MultiHeadAttention, LayerNormalization, Bidirectional
import numpy as np
from configs import config
from configs import gpu_config
import tensorflow as tf

class Hybrid_Model:
    def __init__(self, output_class_num):
        self.class_num = output_class_num

    def train(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def build_model(self):
        # Create a MirroredStrategy.
        strategy = tf.distribute.MirroredStrategy()
        # Open a strategy scope and create/restore the model
        with strategy.scope():
            # Build a hybrid model
            model = self.hybrid_model(self.cnn_model, self.lstm_model)
        return model

    # build CNN model
    def cnn_model(self, input_shape):
        # Input layer
        input_layer = Input(shape=input_shape)
        # Add convolutional layers
        conv1 = Conv1D(32, 3, activation='relu')(input_layer)
        conv2 = Conv1D(32, 3, activation='relu')(conv1)
        conv3 = Conv1D(32, 3, activation='relu')(conv2)
        dropout1 = Dropout(0.5)(conv3)
        # Add flatten and fully connected layers
        flatten = Flatten()(dropout1)
        # Build the model
        model = Model(inputs=input_layer, outputs=flatten)
        return model

    def lstm_model(self, seq_len, num_features, num_lstm_layers, hidden_size, dropout_rate):
        # LSTM model
        input_layer = Input(shape=(seq_len, num_features))
        lstm = input_layer
        for i in range(num_lstm_layers):
            if i == num_lstm_layers - 1:
                lstm = Bidirectional(LSTM(units=hidden_size, dropout=dropout_rate, return_sequences=False))(lstm)
            else:
                lstm = Bidirectional(LSTM(units=hidden_size, dropout=dropout_rate, return_sequences=True))(lstm)
        output_layer = lstm
        model = Model(input_layer, output_layer)
        return model


    # Build the Transformer model
    def transformer_model(self, input_shape, vocab_size, num_heads, d_model, num_layers, max_seq_length):
        inputs = Input(shape=input_shape)
        embedding_layer = Embedding(input_dim=vocab_size, output_dim=d_model)(inputs)

        # Generate positional encodings
        position_encoding = self.positional_encoding(max_seq_length, d_model)
        encoded_inputs = embedding_layer + position_encoding

        # Create multiple Transformer layers
        for _ in range(num_layers):
            attention_output = MultiHeadAttention(num_heads=num_heads, key_dim=d_model)(value=encoded_inputs,
                                                                                        query=encoded_inputs,
                                                                                        key=encoded_inputs)
            attention_output = LayerNormalization(epsilon=1e-6)(attention_output + encoded_inputs)

            feedforward_output = self.feed_forward(attention_output, d_model)
            encoded_inputs = LayerNormalization(epsilon=1e-6)(feedforward_output + attention_output)

        # Perform global average pooling on the final output
        pooling_layer = GlobalAveragePooling1D()(encoded_inputs)

        return Model(inputs=inputs, outputs=pooling_layer)

    # Generate positional encodings
    def positional_encoding(self, max_seq_length, d_model):
        position = np.arange(max_seq_length)[:, np.newaxis]
        angle_rates = 1 / np.power(10000, (2 * (np.arange(d_model)[np.newaxis, :]) // 2) / d_model)
        position_encoding = position * angle_rates

        # Use sine function for even indices and cosine function for odd indices
        position_encoding[:, 0::2] = np.sin(position_encoding[:, 0::2])
        position_encoding[:, 1::2] = np.cos(position_encoding[:, 1::2])

        return position_encoding

    # Feedforward neural network
    def feed_forward(self, x, d_model):
        d_ff = 4 * d_model  # 前馈层维度
        ff_layer = tf.keras.Sequential([
            Dense(d_ff, activation='relu'),
            Dense(d_model)
        ])
        return ff_layer(x)

    # Build and compile a hybrid model
    def hybrid_model(self, cnn_model, lstm_model):
        cnn_input_shape = (config.X_feature_len, 1)
        cnn1 = cnn_model(cnn_input_shape)
        cnn2 = cnn_model(cnn_input_shape)

        # lstm1 = lstm_model(config.lstm_seq_len, config.lstm_features_num, config.lstm_num_layers, config.lstm_hidden_size, config.lstm_dropout)
        # lstm2 = lstm_model(config.lstm_seq_len, config.lstm_features_num, config.lstm_num_layers,
        #                    config.lstm_hidden_size, config.lstm_dropout)

        combined_output = concatenate([cnn1.output, cnn2.output])
        dense_layer = Dense(128, activation='relu')(combined_output)
        output_layer = Dense(int(self.class_num), activation='softmax')(dense_layer)

        model = Model(inputs=[cnn1.input, cnn2.input], outputs=output_layer)
        model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
        return model