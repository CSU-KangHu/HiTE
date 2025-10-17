import torch
import torch.nn as nn
import torch.nn.functional as F

class FocalLoss(nn.Module):
    def __init__(self, alpha=torch.tensor([0.7, 0.3]), gamma=2):
        super(FocalLoss, self).__init__()
        self.alpha = alpha   
        self.gamma = gamma   
        
    def forward(self, inputs, targets):
        self.alpha = self.alpha.to(targets.device)
        ce_loss = F.cross_entropy(inputs, targets, reduction='none')   
        pt = torch.exp(-ce_loss)   
        focal_loss = self.alpha[targets] * (1-pt)**self.gamma * ce_loss 
        return focal_loss.mean()