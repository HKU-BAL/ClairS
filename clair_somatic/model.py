import torch
from torch import nn, einsum
import torch.nn.functional as F

from einops import rearrange, repeat
from einops.layers.torch import Rearrange
import numpy as np
import shared.param as param
# helper methods

def cal_scale(input_size, layers):
    output_size = input_size
    for _ in range(layers):
        output_size = np.ceil(output_size / 2)
    return int(output_size)

import torch
import torch.nn as nn

class BasicBlock(nn.Module):
    """Basic Block for resnet 18 and resnet 34
    """

    #BasicBlock and BottleNeck block
    #have different output size
    #we use class attribute expansion
    #to distinct
    expansion = 1

    def __init__(self, in_channels, out_channels, stride=1):
        super().__init__()

        #residual function
        self.residual_function = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels * BasicBlock.expansion, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(out_channels * BasicBlock.expansion)
        )

        #shortcut
        self.shortcut = nn.Sequential()

        #the shortcut output dimension is not the same with residual function
        #use 1*1 convolution to match the dimension
        if stride != 1 or in_channels != BasicBlock.expansion * out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels * BasicBlock.expansion, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels * BasicBlock.expansion)
            )

    def forward(self, x):
        return nn.ReLU(inplace=True)(self.residual_function(x) + self.shortcut(x))


class ResNet_pileup(nn.Module):

    def __init__(self, apply_softmax=False, block=BasicBlock, num_block=1, num_classes=3, ):
        super().__init__()

        channel_size = param.pileup_channel_size
        self.dim = channel_size * 2
        in_channels = 4
        conv1_channel_size = 16
        conv2_channel_size = 32
        conv3_channel_size = 64

        self.conv1 = BasicConv2D(input_channel=in_channels, output_channel=conv1_channel_size, strides=2)

        depth_scale_size = cal_scale(input_size=self.dim//in_channels, layers=3)
        width_scale_size = cal_scale(input_size=param.no_of_positions, layers=3)

        self.conv1_x = self._make_layer(block, conv1_channel_size, 1, 1)

        self.conv2 = BasicConv2D(input_channel=conv1_channel_size, output_channel=conv2_channel_size, strides=2)
        self.conv2_x = self._make_layer(block, conv2_channel_size, 1, 1)

        self.conv3 = BasicConv2D(input_channel=conv2_channel_size, output_channel=conv3_channel_size, strides=2)
        self.conv3_x = self._make_layer(block, conv3_channel_size, 1, 1)

        self.dropout = nn.Dropout(p=0.5)
        # self.conv3_x = self._make_layer(block, 64, num_block[1], 2)
        # self.conv4_x = self._make_layer(block, 256, num_block[2], 2)
        # self.conv5_x = self._make_layer(block, 512, num_block[3], 2)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(conv3_channel_size * depth_scale_size * width_scale_size, 16)
        self.fc2 = nn.Linear(16, num_classes)

        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)

    def _make_layer(self, block, out_channels, num_blocks, stride):
        """make resnet layers(by layer i didnt mean this 'layer' was the
        same as a neuron netowork layer, ex. conv layer), one layer may
        contain more than one residual block
        Args:
            block: block type, basic block or bottle neck block
            out_channels: output depth channel number of this layer
            num_blocks: how many blocks per layer
            stride: the stride of the first block of this layer
        Return:
            return a resnet layer
        """

        # we have num_block blocks per layer, the first block
        # could be 1 or 2, other blocks would always be 1
        strides = [stride] #+ [1]
        layers = []
        for stride in strides:
            layers.append(block(out_channels, out_channels, stride))
            self.in_channels = out_channels * block.expansion

        return nn.Sequential(*layers)

    def forward(self, x):
        input = torch.reshape(x, (-1, param.no_of_positions, self.dim // 4, 4))
        input = input.permute(0, 3, 1, 2)
        output = self.conv1(input)
        output = self.conv1_x(output)
        output = self.conv2(output)
        output = self.conv2_x(output)
        output = self.conv3(output)
        output = self.conv3_x(output)

        output = self.dropout(output)
        output = self.flatten(output)
        # output = self.conv3_x(output)
        # output = self.conv4_x(output)
        # output = self.conv5_x(output)
        # output = self.avg_pool(output)
        output = output.view(output.size(0), -1)
        output = self.fc1(output)
        output = self.fc2(output)

        if self.apply_softmax:
            output = self.softmax(output)
        return output

class bigru(nn.Module):
    def __init__(
        self,
        num_classes=3,
        width = param.no_of_positions,
        # dim = param.channel_size,
        batch_first=True,
        apply_softmax=False
    ):
        super().__init__()
        kwargs = dict(locals())

        self.num_layers = 2
        self.flatten = nn.Flatten()
        self.lstm_hidden_size = 128
        self.lstm_hidden_size2 = 160
        self.input_shape = [param.no_of_positions, self.lstm_hidden_size2 * 2]
        # self.dropout = nn.Dropout(0.2)

        channel_size = param.pileup_channel_size
        self.dim = channel_size * 2
        self.lstm = nn.LSTM(input_size=self.dim, hidden_size=self.lstm_hidden_size, batch_first=batch_first,
                            num_layers=1, bidirectional=True)

        self.lstm_2 = nn.LSTM(input_size=self.lstm_hidden_size*2, hidden_size=self.lstm_hidden_size2, batch_first=batch_first,
                            num_layers=1, bidirectional=True)

        fc_unit_1 = 128
        self.dropout_fc1 = nn.Dropout(p=0.2)
        # depth_scale_size = cal_scale(input_size=depth, layers=len(self.layers_prefix))
        # width_scale_size = cal_scale(input_size=width, layers=len(self.layers_prefix))
        # fc1_shape = dim * depth_scale_size * width_scale_size
        self.fc1 = nn.Linear(in_features=self.input_shape[0] * self.input_shape[1], out_features=fc_unit_1)
        self.dropout_fc2 = nn.Dropout(p=0.2)

        self.fc2 = nn.Linear(in_features=fc_unit_1, out_features=fc_unit_1)
        self.dropout = nn.Dropout(p=0.5)

        self.fc3 = nn.Linear(in_features=fc_unit_1, out_features=num_classes)
        self.selu = nn.SELU()
        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)
    def forward(self, x):
        # h0 = torch.randn(self.num_layers * 2, self.input_shape[0], self.input_shape[1])
        # c0 = torch.randn(self.num_layers * 2, self.input_shape[0], self.input_shape[1])
        # input = torch.reshape(x,(-1, param.no_of_positions, self.dim // 4, 4))
        output, hidden = self.lstm(x)
        # output = self.flatten(output)
        # output = self.layer2(output)
        # if len(self.layers_prefix) == 3:
        #     output = self.layer3(output)
        output, hidden = self.lstm_2(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        # output = self.dropout(output)
        output = self.dropout_fc1(self.fc1(output))

        output = self.dropout_fc2(self.fc2(output))
        output = self.selu(self.fc3(output))

        if self.apply_softmax:
            output = self.softmax(output)
        return output
        # output = self.fc(output)
        # return output
        # return self.layers(x)


#ResNet in Clair3
import torch
import torch.nn as nn


class BottleNeck(nn.Module):
    """Residual block for resnet over 50 layers
    """
    expansion = 4
    def __init__(self, in_channels, out_channels, stride=1):
        super().__init__()
        self.residual_function = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels, stride=stride, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels * BottleNeck.expansion, kernel_size=1, bias=False),
            nn.BatchNorm2d(out_channels * BottleNeck.expansion),
        )

        self.shortcut = nn.Sequential()

        if stride != 1 or in_channels != out_channels * BottleNeck.expansion:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels * BottleNeck.expansion, stride=stride, kernel_size=1, bias=False),
                nn.BatchNorm2d(out_channels * BottleNeck.expansion)
            )

    def forward(self, x):
        return nn.ReLU(inplace=True)(self.residual_function(x) + self.shortcut(x))



class BasicConv2D(nn.Module):
    def __init__(self, input_channel, output_channel, kernel_size=3, strides=1,padding=1):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv2d(input_channel, output_channel, kernel_size=kernel_size, stride=strides, padding=padding, bias=False),
            nn.BatchNorm2d(output_channel),
            nn.ReLU(inplace=True))

    def forward(self, x):
        return self.conv(x)

class ResNet(nn.Module):

    def __init__(self, block=BasicBlock, num_block=1, num_classes=3):
        super().__init__()

        in_channels = param.channel_size
        conv1_channel_size = 16
        conv2_channel_size = 32
        conv3_channel_size = 64
        self.conv1 = BasicConv2D(input_channel=in_channels, output_channel=conv1_channel_size, strides=2)

        depth_scale_size = cal_scale(input_size=param.max_depth, layers=3)
        width_scale_size = cal_scale(input_size=param.no_of_positions, layers=3)

        self.conv1_x = self._make_layer(block, conv1_channel_size, 1, 1)

        self.conv2 = BasicConv2D(input_channel=conv1_channel_size, output_channel=conv2_channel_size, strides=2)
        self.conv2_x = self._make_layer(block, conv2_channel_size, 1, 1)

        self.conv3 = BasicConv2D(input_channel=conv2_channel_size, output_channel=conv3_channel_size, strides=2)
        self.conv3_x = self._make_layer(block, conv3_channel_size, 1, 1)

        self.dropout = nn.Dropout(p=0.5)
        # self.conv3_x = self._make_layer(block, 64, num_block[1], 2)
        # self.conv4_x = self._make_layer(block, 256, num_block[2], 2)
        # self.conv5_x = self._make_layer(block, 512, num_block[3], 2)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(conv3_channel_size * depth_scale_size * width_scale_size, 16)
        self.fc2 = nn.Linear(16, num_classes)

    def _make_layer(self, block, out_channels, num_blocks, stride):
        """make resnet layers(by layer i didnt mean this 'layer' was the
        same as a neuron netowork layer, ex. conv layer), one layer may
        contain more than one residual block
        Args:
            block: block type, basic block or bottle neck block
            out_channels: output depth channel number of this layer
            num_blocks: how many blocks per layer
            stride: the stride of the first block of this layer
        Return:
            return a resnet layer
        """

        # we have num_block blocks per layer, the first block
        # could be 1 or 2, other blocks would always be 1
        strides = [stride] #+ [1]
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

        output = self.dropout(output)
        output = self.flatten(output)
        # output = self.conv3_x(output)
        # output = self.conv4_x(output)
        # output = self.conv5_x(output)
        # output = self.avg_pool(output)
        output = output.view(output.size(0), -1)
        output = self.fc1(output)
        output = self.fc2(output)

        return output
