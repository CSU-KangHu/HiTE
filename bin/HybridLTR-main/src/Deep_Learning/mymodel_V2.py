import torch
import torch.nn as nn
torch.manual_seed(2024)
import torch.nn.functional as F

class LSTMCat(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(LSTMCat, self).__init__()
        # self.lstm2 = nn.LSTM(1, hidden_dim//2, num_layers=3, batch_first=True, bidirectional=True)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        # 左边
        self.cnn_left = nn.Sequential(
            nn.Conv1d(in_channels=1, out_channels=5, kernel_size=5, stride=1, padding=0),  # 775 - 5 + 1 = 771
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Conv1d(in_channels=5, out_channels=5, kernel_size=5, stride=3, padding=0),  # (771 -5 ) / 3 + 1 = 256
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Conv1d(in_channels=5, out_channels=5, kernel_size=3, stride=1, padding=0), # (256 - 3) / 1 + 1 = 254
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Flatten(start_dim=1)
        )
        # 右边
        self.cnn_right = nn.Sequential(
            nn.Conv1d(in_channels=1, out_channels=5, kernel_size=5, stride=1, padding=0),  # 775 - 5 + 1 = 771
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Conv1d(in_channels=5, out_channels=5, kernel_size=5, stride=3, padding=0),  # (771 -5 ) / 3 + 1 = 256
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Conv1d(in_channels=5, out_channels=5, kernel_size=3, stride=1, padding=0), # (256 - 3) / 1 + 1 = 254
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Flatten(start_dim=1)
        )

        self.seq3 = nn.Sequential(
            nn.Linear(254 * 5 * 2, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.BatchNorm1d(512),
            nn.ReLU()
        )

        self.combine_cnn = nn.Sequential(
            nn.Conv1d(in_channels=2, out_channels=3, kernel_size=3, stride=1, padding=0),  # 512 -3  + 1 = 510
            nn.BatchNorm1d(3),
            nn.ReLU(),
            nn.Conv1d(in_channels=3, out_channels=5, kernel_size=3, stride=1, padding=0),  # 510-3 + 1 = 508
            nn.BatchNorm1d(5),
            nn.ReLU(),
            nn.Conv1d(in_channels=5, out_channels=3, kernel_size=3, stride=1, padding=0), # 508 - 3 + 1 = 506
            nn.BatchNorm1d(3),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Flatten(start_dim=1)
        )

        self.seq4 = nn.Sequential(
            nn.Linear(506 * 3, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Linear(64, output_dim)
        )

    def forward(self, x, x2):
        lx = x2[:, :1, :]
        rx = x2[:, 1:, :]
        cnn_out1 = self.cnn_left(lx)
        cnn_out2 = self.cnn_right(rx)
        cnn_out = torch.cat((cnn_out1, cnn_out2), 1)
        # print(cnn_out.shape)
        cnn_out = self.seq3(cnn_out)
        # print(cnn_out.shape)
        x = x.permute(0, 2, 1)
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # print(output.shape)
        output = torch.stack([output, cnn_out], dim=1)
        # print(output.shape)
        output = self.combine_cnn(output)
        output = self.seq4(output)

        return output