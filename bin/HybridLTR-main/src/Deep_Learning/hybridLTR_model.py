import torch
import torch.nn as nn
torch.manual_seed(2233)
import torch.nn.functional as F

class GlobalAvgPool2d(nn.Module):
    # define global average pooling layer
    def __init__(self):
        super(GlobalAvgPool2d, self).__init__()
    
    def forward(self, x):
        return F.adaptive_avg_pool2d(x, (1, 1))

class BasicBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, dropout_rate=0.2):
        super(BasicBlock, self).__init__()
        self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.dropout1 = nn.Dropout2d(p=dropout_rate)
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)
        self.dropout2 = nn.Dropout2d(p=dropout_rate)

        self.shortcut = nn.Sequential()
        if stride != 1 or in_channels != out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels)
            )

    def forward(self, x):
        out = self.conv1(x)
        out = self.bn1(out)
        out = F.relu(out)
        out = self.dropout1(out)
        
        out = self.conv2(out)
        out = self.bn2(out)
        out = self.dropout2(out)
        
        out += self.shortcut(x)
        out = F.relu(out)
        return out

    
class CNNCAT(nn.Module): # (batch_size, channels, height, width)
    def __init__(self, dropout_rate=0.3):
        super(CNNCAT, self).__init__()
        self.model = nn.Sequential(
            BasicBlock(3, 8, dropout_rate=dropout_rate),
            BasicBlock(8, 16, dropout_rate=dropout_rate),
            BasicBlock(16, 32, dropout_rate=dropout_rate)
        )

        self.model_kmer = nn.Sequential(
            BasicBlock(2, 8, dropout_rate=dropout_rate),
            BasicBlock(8, 16, dropout_rate=dropout_rate),
            BasicBlock(16, 32, dropout_rate=dropout_rate)
        )

        self.global_avg_pool = GlobalAvgPool2d()

        self.flatten = nn.Flatten()

        self.fc_layers = nn.Sequential(
            nn.Linear(32 * 2, 16),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Dropout(dropout_rate),

            nn.Linear(16, 8),
            nn.BatchNorm1d(8),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
        
            nn.Linear(8, 2)
        )

    def forward(self, x, x_kmer):
        out = self.model(x)
        out_kmer = self.model_kmer(x_kmer)

        out = self.global_avg_pool(out)
        out_kmer = self.global_avg_pool(out_kmer)

        out = self.flatten(out)
        out_kmer = self.flatten(out_kmer)

        out = torch.cat((out, out_kmer), 1)
        out = self.fc_layers(out)
        
        return out