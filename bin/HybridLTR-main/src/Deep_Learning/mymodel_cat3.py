import torch
import torch.nn as nn
torch.manual_seed(2233)
import torch.nn.functional as F

class GlobalAvgPool2d(nn.Module):
    # 定义全局平均池化层
    def __init__(self):
        super(GlobalAvgPool2d, self).__init__()
    
    def forward(self, x):
        return F.adaptive_avg_pool2d(x, (1, 1))

class BasicBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1):
        super(BasicBlock, self).__init__()
        self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)

        self.shortcut = nn.Sequential()
        if stride != 1 or in_channels != out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels)
            )

    def forward(self, x):
        out = torch.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        out += self.shortcut(x)
        out = torch.relu(out)
        return out

def make_conv_block(in_channels, out_channels, kernel_size=3, padding=1):
    return nn.Sequential(
        nn.Conv2d(in_channels=in_channels, out_channels=out_channels, kernel_size=kernel_size, padding=padding),
        nn.BatchNorm2d(out_channels),
        nn.ReLU()
    )
    
# 输入维度：[N，4, 100, 200]
class LSTMCat(nn.Module): # (batch_size, channels, height, width)
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(LSTMCat, self).__init__()
        self.model = nn.Sequential(
            BasicBlock(3, 8),
            BasicBlock(8, 32),
            BasicBlock(32, 64)
        )

        self.model_kmer = nn.Sequential(
            BasicBlock(2, 8),
            BasicBlock(8, 32),
            BasicBlock(32, 64)
        )

        self.global_avg_pool = GlobalAvgPool2d()

        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(64 * 2, 32)
        self.fc2 = nn.Linear(32, 2)

    def forward(self, x, x_kmer):
        out = self.model(x)
        out_kmer = self.model_kmer(x_kmer)

        out = self.global_avg_pool(out)
        out_kmer = self.global_avg_pool(out_kmer)

        out = self.flatten(out)
        out_kmer = self.flatten(out_kmer)

        out = torch.cat((out, out_kmer), 1)

        out = self.fc1(out)
        out = self.fc2(out)
        
        return out