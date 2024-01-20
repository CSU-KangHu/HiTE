import os

from keras.models import load_model, Model
from keras.layers import Input, Dense, Dropout, Flatten, Conv1D,\
    Embedding, GlobalAveragePooling1D, concatenate, \
    MultiHeadAttention, LayerNormalization
import numpy as np
from configs import config
from configs import gpu_config
import tensorflow as tf

class Hybrid_Model:
    def __init__(self, input_features_num, output_class_num):
        self.num_features = input_features_num
        self.class_num = output_class_num

    def train(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def build_model(self):
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

            # Build a hybrid model
            model = self.hybrid_model(self.cnn_model, self.transformer_model)
        return model

    # build CNN model
    def cnn_model(self, input_shape):
        # Input layer
        input_layer = Input(shape=input_shape)
        # Add convolutional layers
        conv1 = Conv1D(32, 3, activation='relu')(input_layer)
        conv2 = Conv1D(32, 3, activation='relu')(conv1)
        conv3 = Conv1D(32, 3, activation='relu')(conv2)
        conv4 = Conv1D(32, 3, activation='relu')(conv3)
        dropout1 = Dropout(0.5)(conv4)
        # Add flatten and fully connected layers
        flatten = Flatten()(dropout1)
        # Build the model
        model = Model(inputs=input_layer, outputs=flatten)
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
    def hybrid_model(self, cnn_model, transformer_model):
        cnn_input_shape = (self.num_features, 1)
        max_TSD_seq_length = 11
        vocab_size = 5
        num_heads = 4
        d_model = 16
        num_layers = 2
        transformer_input_shape = (max_TSD_seq_length,)

        cnn = cnn_model(cnn_input_shape)
        transformer = transformer_model(transformer_input_shape, vocab_size, num_heads, d_model, num_layers,
                                        max_TSD_seq_length)

        combined_output = concatenate([cnn.output, transformer.output])
        dense_layer = Dense(128, activation='relu')(combined_output)
        # dense_layer = Dense(128, activation='relu')(cnn.output)
        output_layer = Dense(int(self.class_num), activation='softmax')(dense_layer)

        model = Model(inputs=[cnn.input, transformer.input], outputs=output_layer)
        # model = Model(inputs=cnn.input, outputs=output_layer)
        model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
        return model