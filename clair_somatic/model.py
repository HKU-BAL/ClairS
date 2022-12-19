import numpy as np
from torch import nn

import shared.param as param


def cal_scale(input_size, layers):
    output_size = input_size
    for _ in range(layers):
        output_size = np.ceil(output_size / 2)
    return int(output_size)


class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, in_channels, out_channels, stride=1):
        super().__init__()

        # residual function
        self.residual_function = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels * BasicBlock.expansion, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(out_channels * BasicBlock.expansion)
        )

        # shortcut
        self.shortcut = nn.Sequential()

        if stride != 1 or in_channels != BasicBlock.expansion * out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels * BasicBlock.expansion, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels * BasicBlock.expansion)
            )

    def forward(self, x):
        return nn.ReLU(inplace=True)(self.residual_function(x) + self.shortcut(x))


class bigru(nn.Module):
    def __init__(
            self,
            num_classes=3,
            width=param.no_of_positions,
            batch_first=True,
            apply_softmax=False,
            channel_size=param.pileup_channel_size * 2
    ):
        super().__init__()
        kwargs = dict(locals())

        self.num_layers = 2
        self.flatten = nn.Flatten()
        self.lstm_hidden_size = 128
        self.lstm_hidden_size2 = 192
        fc1_layer_size = 128
        dropout_rate = 0.3

        self.input_shape = [param.no_of_positions, self.lstm_hidden_size2 * 2]

        self.dim = channel_size
        self.lstm = nn.GRU(input_size=self.dim, hidden_size=self.lstm_hidden_size, batch_first=batch_first,
                           num_layers=1, bidirectional=True)

        self.lstm_2 = nn.GRU(input_size=self.lstm_hidden_size * 2, hidden_size=self.lstm_hidden_size2,
                             batch_first=batch_first,
                             num_layers=1, bidirectional=True)

        self.dropout_fc1 = nn.Dropout(p=dropout_rate)
        self.fc1 = nn.Linear(in_features=self.input_shape[0] * self.input_shape[1], out_features=fc1_layer_size)
        self.dropout_fc2 = nn.Dropout(p=dropout_rate)

        self.fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)

        self.fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.selu = nn.SELU()
        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)

    def forward(self, x):

        output, hidden = self.lstm(x)
        output, hidden = self.lstm_2(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout_fc1(output)
        output = self.selu(self.dropout_fc1(self.fc1(output)))

        output = self.selu(self.dropout_fc2(self.fc2(output)))
        output = self.selu(self.fc3(output))

        if self.apply_softmax:
            output = self.softmax(output)
        return output


class BasicConv2D(nn.Module):
    def __init__(self, input_channel, output_channel, kernel_size=3, strides=1, padding=1):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv2d(input_channel, output_channel, kernel_size=kernel_size, stride=strides, padding=padding,
                      bias=False),
            nn.BatchNorm2d(output_channel),
            nn.ReLU(inplace=True))

    def forward(self, x):
        return self.conv(x)


class ResNet(nn.Module):

    def __init__(self, block=BasicBlock, num_block=1, platform='ont', num_classes=3):
        super().__init__()

        in_channels = param.channel_size
        conv1_channel_size = 64
        conv2_channel_size = 96
        conv3_channel_size = 128
        fc1_layer_size = 160
        fc2_layer_size = 128
        dropout_rate = 0.3

        self.conv1 = BasicConv2D(input_channel=in_channels, output_channel=conv1_channel_size, strides=2)

        depth_scale_size = cal_scale(input_size=param.matrix_depth_dict[platform], layers=3)
        width_scale_size = cal_scale(input_size=param.no_of_positions, layers=3)

        self.conv1_x = self._make_layer(block, conv1_channel_size)

        self.conv2 = BasicConv2D(input_channel=conv1_channel_size, output_channel=conv2_channel_size, strides=2)
        self.conv2_x = self._make_layer(block, conv2_channel_size)

        self.conv3 = BasicConv2D(input_channel=conv2_channel_size, output_channel=conv3_channel_size, strides=2)
        self.conv3_x = self._make_layer(block, conv3_channel_size)

        self.dropout = nn.Dropout(p=dropout_rate)
        self.selu = nn.SELU()

        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(conv3_channel_size * depth_scale_size * width_scale_size, fc1_layer_size)
        self.fc2 = nn.Linear(fc1_layer_size, fc2_layer_size)
        self.fc3 = nn.Linear(fc2_layer_size, num_classes)

    def _make_layer(self, block, out_channels, stride=1):
        strides = [stride]
        layers = []
        for stride in strides:
            layers.append(block(out_channels, out_channels, stride))
            self.in_channels = out_channels * block.expansion

        return nn.Sequential(*layers)

    def forward(self, x):
        output = self.conv1(x)
        output = self.conv1_x(output)
        output = self.conv2(output)
        output = self.conv2_x(output)
        output = self.conv3(output)
        output = self.conv3_x(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout(output)
        output = self.selu(self.dropout(self.fc1(output)))
        output = self.selu(self.dropout(self.fc2(output)))
        output = self.selu(self.fc3(output))

        return output
