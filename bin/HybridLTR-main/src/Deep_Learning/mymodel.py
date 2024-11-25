import torch
import torch.nn as nn
torch.manual_seed(2024)
import torch.nn.functional as F

class AttentionLSTM(nn.Module):
    def __init__(self, input_dim=9, hidden_dim=512, output_dim=2):
        super(AttentionLSTM, self).__init__()
        # self.embedding = nn.Embedding(100, input_dim)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        # nn.init.xavier_normal_(self.lstm.weight_ih_l0)
        self.fc1 = nn.Linear(hidden_dim, 32)
        nn.init.xavier_normal_(self.fc1.weight)
        self.fc2 = nn.Linear(32, output_dim)
        nn.init.xavier_normal_(self.fc2.weight)
        # self.softmax = nn.Softmax(dim=1)

    def attention(self, lstm_output, final_state):
        # lstm_output = lstm_output.permute(1, 0, 2)
        # print(lstm_output.shape) #([batch_size, seq_length, hidden_dim*2])
        # merged_state = torch.cat([s for s in final_state], 1)
        merged_state = final_state
        # print(merged_state.shape) # torch.Size([batch_size, hidden_dim * num_layers])
        merged_state = merged_state.squeeze(0).unsqueeze(2)
        # print(merged_state.shape) # torch.Size([batch_size, hidden_dim * num_layers, 1])
        weights = torch.bmm(lstm_output, merged_state)
        weights = F.softmax(weights.squeeze(2), dim=1).unsqueeze(2)
        return torch.bmm(torch.transpose(lstm_output, 1, 2), weights).squeeze(2)

    def forward(self, x):
        # embedded = self.embedding(x)
        # print(embedded.shape)
        output, (hidden, cell) = self.lstm(x)
        # num_layers, num_directions, batch, hidden_size
        hidden = hidden.view(5, 2, -1, 512//2)
        final_hidden = hidden[-1]
        final_hidden = final_hidden.permute(1, 0, 2)
        final_hidden = final_hidden.reshape(final_hidden.shape[0], -1)
        # print(final_hidden.shape)        
        # new_hidden = hidden.view(num_layers, num_directions, batch, hidden_size)
        # attention_weights = self.softmax(output[:,-1,:])
        # attended_output = torch.sum(attention_weights * output, dim=1)
        # output = self.fc1(attended_output)
        attn_output = self.attention(output, final_hidden)
        # print(attn_output.shape) #torch.Size([32, 128])
        output = self.fc1(attn_output.squeeze(0)) 
        # print(output.shape) # ([32, 100, 32])
        output = F.relu(output)
        output = self.fc2(output)
        # print(output.shape) # torch.Size([32, 100, 2])
        return output

class AttentionLSTM5(nn.Module):
    def __init__(self, input_dim=5, hidden_dim=512, output_dim=2):
        super(AttentionLSTM5, self).__init__()
        # self.embedding = nn.Embedding(100, input_dim)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        # nn.init.xavier_normal_(self.lstm.weight_ih_l0)
        self.fc1 = nn.Linear(hidden_dim, 32)
        nn.init.xavier_normal_(self.fc1.weight)
        self.fc2 = nn.Linear(32, output_dim)
        nn.init.xavier_normal_(self.fc2.weight)
        # self.softmax = nn.Softmax(dim=1)

    def attention(self, lstm_output, final_state):
        # lstm_output = lstm_output.permute(1, 0, 2)
        # print(lstm_output.shape) #([batch_size, seq_length, hidden_dim*2])
        # merged_state = torch.cat([s for s in final_state], 1)
        merged_state = final_state
        # print(merged_state.shape) # torch.Size([batch_size, hidden_dim * num_layers])
        merged_state = merged_state.squeeze(0).unsqueeze(2)
        # print(merged_state.shape) # torch.Size([batch_size, hidden_dim * num_layers, 1])
        weights = torch.bmm(lstm_output, merged_state)
        weights = F.softmax(weights.squeeze(2), dim=1).unsqueeze(2)
        return torch.bmm(torch.transpose(lstm_output, 1, 2), weights).squeeze(2)

    def forward(self, x):
        # embedded = self.embedding(x)
        # print(embedded.shape)
        output, (hidden, cell) = self.lstm(x)
        hidden = hidden.view(5, 2, -1, 512//2)
        final_hidden = hidden[-1]
        final_hidden = final_hidden.permute(1, 0, 2)
        final_hidden = final_hidden.reshape(final_hidden.shape[0], -1)
        # print(final_hidden.shape)        
        # new_hidden = hidden.view(num_layers, num_directions, batch, hidden_size)
        # attention_weights = self.softmax(output[:,-1,:])
        # attended_output = torch.sum(attention_weights * output, dim=1)
        # output = self.fc1(attended_output)
        attn_output = self.attention(output, final_hidden)
        # print(attn_output.shape) #torch.Size([32, 128])
        output = self.fc1(attn_output.squeeze(0)) 
        # print(output.shape) # ([32, 100, 32])
        output = F.relu(output)
        output = self.fc2(output)
        # print(output.sha
        return output

class AttentionLSTM8(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTM8, self).__init__()
        self.input_normalization = nn.BatchNorm1d(100)
        # self.embedding = nn.Embedding(100, input_dim)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        # nn.init.xavier_normal_(self.lstm.weight_ih_l0)
        self.dropout = nn.Dropout(0.1)

        self.fc1 = nn.Linear(hidden_dim, 32)
        nn.init.xavier_normal_(self.fc1.weight)
        self.fc2 = nn.Linear(32, output_dim)
        nn.init.xavier_normal_(self.fc2.weight)
        # self.softmax = nn.Softmax(dim=1)

    def attention(self, lstm_output, final_state):
        # lstm_output = lstm_output.permute(1, 0, 2)
        # print(lstm_output.shape) #([batch_size, seq_length, hidden_dim*2])
        # merged_state = torch.cat([s for s in final_state], 1)
        merged_state = final_state
        # print(merged_state.shape) # torch.Size([batch_size, hidden_dim * num_layers])
        merged_state = merged_state.squeeze(0).unsqueeze(2)
        # print(merged_state.shape) # torch.Size([batch_size, hidden_dim * num_layers, 1])
        weights = torch.bmm(lstm_output, merged_state)
        weights = F.softmax(weights.squeeze(2), dim=1).unsqueeze(2)
        return torch.bmm(torch.transpose(lstm_output, 1, 2), weights).squeeze(2)

    def forward(self, x):
        # embedded = self.embedding(x)
        # print(embedded.shape)
        x = self.input_normalization(x)
        output, (hidden, cell) = self.lstm(x)
        output = self.dropout(output)
        hidden = hidden.view(5, 2, -1, 512//2)
        final_hidden = hidden[-1]
        final_hidden = final_hidden.permute(1, 0, 2)
        final_hidden = final_hidden.reshape(final_hidden.shape[0], -1)
        # print(final_hidden.shape)        
        # new_hidden = hidden.view(num_layers, num_directions, batch, hidden_size)
        # attention_weights = self.softmax(output[:,-1,:])
        # attended_output = torch.sum(attention_weights * output, dim=1)
        # output = self.fc1(attended_output)
        attn_output = self.attention(output, final_hidden)
        # print(attn_output.shape) #torch.Size([32, 128])
        output = self.fc1(attn_output.squeeze(0)) 
        # print(output.shape) # ([32, 100, 32])
        output = F.relu(output)
        output = self.fc2(output)
        # print(output.shape) # torch.Size([32, 100, 2])
        return output
    
class AttentionLSTMCNN(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTMCNN, self).__init__()
        self.input_normalization = nn.BatchNorm1d(8)
        self.conv1 = nn.Conv1d(in_channels=input_dim, out_channels=8, kernel_size=3, stride=1, padding=0)  # (100 - 3 + 2 * 0 ) / 1 + 1 = 98
        nn.init.xavier_normal_(self.conv1.weight)
        self.conv2 = nn.Conv1d(in_channels=8, out_channels=32, kernel_size=3, stride=1, padding=0)  # (98 - 3 + 2 * 0 ) / 1 + 1 = 96
        nn.init.xavier_normal_(self.conv2.weight)
        self.conv3 = nn.Conv1d(in_channels=32, out_channels=8, kernel_size=3, stride=1, padding=0)  # (96 - 3 + 2 * 0 ) / 1 + 1 = 94 [N, 8, 94]
        nn.init.xavier_normal_(self.conv3.weight)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        self.dropout = nn.Dropout(0.2)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight)
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight)
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight)
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)

    def attention(self, lstm_output, final_state):
        merged_state = final_state
        merged_state = merged_state.squeeze(0).unsqueeze(2)
        weights = torch.bmm(lstm_output, merged_state)
        weights = F.softmax(weights.squeeze(2), dim=1).unsqueeze(2)
        return torch.bmm(torch.transpose(lstm_output, 1, 2), weights).squeeze(2)

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = F.relu(self.bn1(x))
        x = self.conv2(x)
        x = F.relu(self.bn2(x))
        x = self.conv3(x)
        x = F.relu(self.bn3(x))
        x = x.permute(0, 2, 1)
        output, (hidden, cell) = self.lstm(x)
        output = self.dropout(output)
        hidden = hidden.view(5, 2, -1, 512//2)
        final_hidden = hidden[-1]
        final_hidden = final_hidden.permute(1, 0, 2)
        final_hidden = final_hidden.reshape(final_hidden.shape[0], -1)
        attn_output = self.attention(output, final_hidden)
        output = self.fc1(attn_output.squeeze(0)) 
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        return output

class AttentionLSTMCNN_BCE(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=1):
        super(AttentionLSTMCNN_BCE, self).__init__()
        self.input_normalization = nn.BatchNorm1d(8)
        self.conv1 = nn.Conv1d(in_channels=input_dim, out_channels=8, kernel_size=3, stride=1, padding=1)  # (n - k + 2 * p ) / s + 1 # (100 - 3 + 2 * 1 ) / 1 + 1 = 98
        nn.init.xavier_normal_(self.conv1.weight)
        self.conv2 = nn.Conv1d(in_channels=8, out_channels=32, kernel_size=3, stride=1, padding=1)  # (98 - 3 + 2 * 1 ) / 1 + 1 = 96
        nn.init.xavier_normal_(self.conv2.weight)
        self.conv3 = nn.Conv1d(in_channels=32, out_channels=8, kernel_size=3, stride=1, padding=1)  # (96 - 3 + 2 * 1 ) / 1 + 1 = 94 [N, 8, 94]
        nn.init.xavier_normal_(self.conv3.weight)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        # self.dropout = nn.Dropout(0.1)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight, gain=nn.init.calculate_gain('relu'))
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight, gain=nn.init.calculate_gain('relu'))
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)
        self.flat = nn.Flatten() # pytorch不会隐式地调整输入的形状

    def attention(self, lstm_output, final_state):
        merged_state = final_state
        merged_state = merged_state.squeeze(0).unsqueeze(2)
        weights = torch.bmm(lstm_output, merged_state)
        weights = F.softmax(weights.squeeze(2), dim=1).unsqueeze(2)
        return torch.bmm(torch.transpose(lstm_output, 1, 2), weights).squeeze(2)

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = F.relu(self.bn1(x))
        x = self.conv2(x)
        x = F.relu(self.bn2(x))
        x = self.conv3(x)
        x = F.relu(self.bn3(x))
        x = x.permute(0, 2, 1)
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # output = self.dropout(output[:, -1, :])
        # hidden = hidden.view(5, 2, -1, 512//2)
        # final_hidden = hidden[-1]
        # final_hidden = final_hidden.permute(1, 0, 2)
        # final_hidden = final_hidden.reshape(final_hidden.shape[0], -1)
        # attn_output = self.attention(output, final_hidden)
        # output = self.fc1(attn_output.squeeze(0)) 
        # output = self.flat(output)
        output = self.fc1(output)
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        # output = self.flat(output)
        # output = output.squeeze()
        # output = torch.nn.softmax(output, dim=1)
        return output
    
class AttentionLSTMCNN(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTMCNN, self).__init__()
        self.input_normalization = nn.BatchNorm1d(8)
        self.conv1 = nn.Conv1d(in_channels=input_dim, out_channels=8, kernel_size=3, stride=1, padding=1)  # (n - k + 2 * p ) / s + 1 # (100 - 3 + 2 * 1 ) / 1 + 1 = 98
        nn.init.xavier_normal_(self.conv1.weight)
        self.conv2 = nn.Conv1d(in_channels=8, out_channels=32, kernel_size=3, stride=1, padding=1)  # (98 - 3 + 2 * 1 ) / 1 + 1 = 96
        nn.init.xavier_normal_(self.conv2.weight)
        self.conv3 = nn.Conv1d(in_channels=32, out_channels=8, kernel_size=3, stride=1, padding=1)  # (96 - 3 + 2 * 1 ) / 1 + 1 = 94 [N, 8, 94]
        nn.init.xavier_normal_(self.conv3.weight)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=5, batch_first=True, bidirectional=True)
        # self.dropout = nn.Dropout(0.1)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight, gain=nn.init.calculate_gain('relu'))
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight, gain=nn.init.calculate_gain('relu'))
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)
        self.flat = nn.Flatten() # pytorch不会隐式地调整输入的形状

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = F.relu(self.bn1(x))
        x = self.conv2(x)
        x = F.relu(self.bn2(x))
        x = self.conv3(x)
        x = F.relu(self.bn3(x))
        x = x.permute(0, 2, 1)
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # output = self.dropout(output[:, -1, :])
        output = self.fc1(output)
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        return output
    
class AttentionLSTMCNN2(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTMCNN2, self).__init__()
        # self.input_normalization = nn.BatchNorm1d(8)
        self.conv1 = nn.Conv1d(in_channels=input_dim, out_channels=8, kernel_size=3, stride=1, padding=1)  # (n - k + 2 * p ) / s + 1 # (100 - 3 + 2 * 1 ) / 1 + 1 = 98
        nn.init.xavier_normal_(self.conv1.weight)
        self.conv2 = nn.Conv1d(in_channels=8, out_channels=32, kernel_size=3, stride=1, padding=1)  # (98 - 3 + 2 * 1 ) / 1 + 1 = 96
        nn.init.xavier_normal_(self.conv2.weight)
        self.conv3 = nn.Conv1d(in_channels=32, out_channels=8, kernel_size=3, stride=1, padding=1)  # (96 - 3 + 2 * 1 ) / 1 + 1 = 94 [N, 8, 94]
        nn.init.xavier_normal_(self.conv3.weight)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=2, batch_first=True, bidirectional=True)
        self.dropout = nn.Dropout(0.2)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight, gain=nn.init.calculate_gain('relu'))
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight, gain=nn.init.calculate_gain('relu'))
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)
        self.flat = nn.Flatten() # pytorch不会隐式地调整输入的形状

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = F.relu(self.bn1(x))
        x = self.conv2(x)
        x = F.relu(self.bn2(x))
        x = self.conv3(x)
        x = F.relu(self.bn3(x))
        x = x.permute(0, 2, 1)
        x = self.dropout(x) # 输入之前进行 dropout
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # output = self.dropout(output[:, -1, :])
        output = self.fc1(output)
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        return output


class AttentionLSTMCNN3(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTMCNN3, self).__init__()
        # self.input_normalization = nn.BatchNorm1d(8)
        self.conv1 = nn.Conv1d(in_channels=input_dim, out_channels=8, kernel_size=3, stride=1, padding=1)  # (n - k + 2 * p ) / s + 1 # (100 - 3 + 2 * 1 ) / 1 + 1 = 98
        nn.init.xavier_normal_(self.conv1.weight)
        self.conv2 = nn.Conv1d(in_channels=8, out_channels=32, kernel_size=3, stride=1, padding=1)  # (98 - 3 + 2 * 1 ) / 1 + 1 = 96
        nn.init.xavier_normal_(self.conv2.weight)
        self.conv3 = nn.Conv1d(in_channels=32, out_channels=8, kernel_size=3, stride=1, padding=1)  # (96 - 3 + 2 * 1 ) / 1 + 1 = 94 [N, 8, 94]
        nn.init.xavier_normal_(self.conv3.weight)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=3, batch_first=True, bidirectional=True)
        self.dropout = nn.Dropout(0.2)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight, gain=nn.init.calculate_gain('relu'))
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight, gain=nn.init.calculate_gain('relu'))
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)
        self.flat = nn.Flatten() # pytorch不会隐式地调整输入的形状

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = F.relu(self.bn1(x))
        x = self.conv2(x)
        x = F.relu(self.bn2(x))
        x = self.conv3(x)
        x = F.relu(self.bn3(x))
        x = x.permute(0, 2, 1)
        x = self.dropout(x) # 输入之前进行 dropout
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # output = self.dropout(output[:, -1, :])
        output = self.fc1(output)
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        return output
    
class AttentionLSTMno(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTMno, self).__init__()
        # self.input_normalization = nn.BatchNorm1d(8)
       
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=3, batch_first=True, bidirectional=True)
        self.dropout = nn.Dropout(0.2)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight, gain=nn.init.calculate_gain('relu'))
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight, gain=nn.init.calculate_gain('relu'))
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)
        self.flat = nn.Flatten() # pytorch不会隐式地调整输入的形状

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
       
        x = x.permute(0, 2, 1)
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # output = self.dropout(output[:, -1, :])
        output = self.fc1(output)
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        return output
    
class AttentionLSTMCNN4(nn.Module):
    def __init__(self, input_dim=8, hidden_dim=512, output_dim=2):
        super(AttentionLSTMCNN4, self).__init__()
        # self.input_normalization = nn.BatchNorm1d(8)
        self.conv1 = nn.Conv1d(in_channels=input_dim, out_channels=8, kernel_size=3, stride=1, padding=1)  # (n - k + 2 * p ) / s + 1 # (100 - 3 + 2 * 1 ) / 1 + 1 = 98
        nn.init.xavier_normal_(self.conv1.weight)
        self.conv2 = nn.Conv1d(in_channels=8, out_channels=32, kernel_size=3, stride=1, padding=1)  # (98 - 3 + 2 * 1 ) / 1 + 1 = 96
        nn.init.xavier_normal_(self.conv2.weight)
        self.conv3 = nn.Conv1d(in_channels=32, out_channels=8, kernel_size=3, stride=1, padding=1)  # (96 - 3 + 2 * 1 ) / 1 + 1 = 94 [N, 8, 94]
        nn.init.xavier_normal_(self.conv3.weight)
        self.lstm = nn.LSTM(input_dim, hidden_dim//2, num_layers=3, batch_first=True, bidirectional=True)
        self.dropout = nn.Dropout(0.2)
        self.fc1 = nn.Linear(hidden_dim, 64)
        nn.init.xavier_normal_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))
        self.fc2 = nn.Linear(64, 16)
        nn.init.xavier_normal_(self.fc2.weight, gain=nn.init.calculate_gain('relu'))
        self.fc3 = nn.Linear(16, output_dim)
        nn.init.xavier_normal_(self.fc3.weight, gain=nn.init.calculate_gain('relu'))
        self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn2 = nn.BatchNorm1d(32)
        self.bn3 = nn.BatchNorm1d(input_dim)
        self.flat = nn.Flatten() # pytorch不会隐式地调整输入的形状

    def forward(self, x):
        # x = self.input_normalization(x)
        # x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = F.relu(self.bn1(x))
        x = self.conv2(x)
        x = F.relu(self.bn2(x))
        x = self.conv3(x)
        x = F.relu(self.bn3(x))
        x = x.permute(0, 2, 1)
        x = self.dropout(x) # 输入之前进行 dropout
        output, (hidden, cell) = self.lstm(x)
        output = output[:, -1, :]
        # output = self.dropout(output[:, -1, :])
        output = self.fc1(output)
        output = F.relu(output)
        output = self.fc2(output)
        output = F.relu(output)
        output = self.fc3(output)
        return output