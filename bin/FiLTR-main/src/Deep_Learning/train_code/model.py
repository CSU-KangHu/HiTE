import torch
import torch.nn as nn
import torch.nn.functional as F

# Set a seed for reproducibility of random initializations
torch.manual_seed(2233)

class GlobalAvgPool2d(nn.Module):
    """
    Defines a global average pooling layer.
    This layer averages the feature maps across the spatial dimensions (height and width),
    resulting in a single value per feature map.
    """
    def __init__(self):
        super(GlobalAvgPool2d, self).__init__()
    
    def forward(self, x):
        return F.adaptive_avg_pool2d(x, (1, 1))

class BasicBlock(nn.Module):
    """
    A basic residual block for the CNN.
    It consists of two convolutional layers with batch normalization, ReLU activation,
    and dropout, along with a shortcut connection.
    """
    def __init__(self, in_channels, out_channels, stride=1, dropout_rate=0.2):
        super(BasicBlock, self).__init__()
        # First convolutional layer
        self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.dropout1 = nn.Dropout2d(p=dropout_rate)
        
        # Second convolutional layer
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)
        self.dropout2 = nn.Dropout2d(p=dropout_rate)

        # Shortcut connection to match dimensions if necessary
        self.shortcut = nn.Sequential()
        if stride != 1 or in_channels != out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels)
            )

    def forward(self, x):
        # Main path
        out = self.conv1(x)
        out = self.bn1(out)
        out = F.relu(out)
        out = self.dropout1(out)
        
        out = self.conv2(out)
        out = self.bn2(out)
        out = self.dropout2(out)
        
        # Add the shortcut connection
        out += self.shortcut(x)
        out = F.relu(out)
        return out

class CNNCAT(nn.Module):
    """
    The main CNN model which processes two types of input features (e.g., image-like and k-mer frequency)
    and concatenates them before feeding them to fully connected layers for classification.
    
    Expected input shape for x: (batch_size, 3, height, width)
    Expected input shape for x_kmer: (batch_size, 2, height, width)
    """
    def __init__(self, dropout_rate=0.3):
        super(CNNCAT, self).__init__()
        # CNN backbone for the first input type (e.g., image features)
        self.model = nn.Sequential(
            BasicBlock(3, 8, dropout_rate=dropout_rate),
            BasicBlock(8, 16, dropout_rate=dropout_rate),
            BasicBlock(16, 32, dropout_rate=dropout_rate)
        )

        # CNN backbone for the second input type (e.g., k-mer features)
        self.model_kmer = nn.Sequential(
            BasicBlock(2, 8, dropout_rate=dropout_rate),
            BasicBlock(8, 16, dropout_rate=dropout_rate),
            BasicBlock(16, 32, dropout_rate=dropout_rate)
        )

        # Pooling and flattening layers
        self.global_avg_pool = GlobalAvgPool2d()
        self.flatten = nn.Flatten()

        # Fully connected layers for classification
        self.fc_layers = nn.Sequential(
            nn.Linear(32 * 2, 16), # 32 from each branch
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Dropout(dropout_rate),

            nn.Linear(16, 8),
            nn.BatchNorm1d(8),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
        
            nn.Linear(8, 2) # Output logits for 2 classes
        )

    def forward(self, x, x_kmer):
        # Manually calculate L2 regularization for all model parameters
        l2_reg = torch.tensor(0., requires_grad=True, device=x.device)
        for param in self.parameters():
            l2_reg = l2_reg + torch.norm(param)
        l2_reg = l2_reg.unsqueeze(0) # Convert scalar to a single-element tensor

        # Process both inputs through their respective CNN backbones
        out = self.model(x)
        out_kmer = self.model_kmer(x_kmer)

        # Apply global average pooling
        out = self.global_avg_pool(out)
        out_kmer = self.global_avg_pool(out_kmer)

        # Flatten the pooled features
        out = self.flatten(out)
        out_kmer = self.flatten(out_kmer)

        # Concatenate the features from both branches
        out = torch.cat((out, out_kmer), 1)
        
        # Pass through the fully connected layers
        out = self.fc_layers(out)
        
        return out, l2_reg