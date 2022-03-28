import torch
from torch import nn, einsum
import torch.nn.functional as F

from einops import rearrange, repeat
from einops.layers.torch import Rearrange
import numpy as np
import shared.param as param
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

        self.conv1 = BasicConv2D(input_channel=in_channels,output_channel=64, strides=2)
        self.in_channels = 32

        depth_scale_size = cal_scale(input_size=param.max_depth, layers=3)
        width_scale_size = cal_scale(input_size=param.no_of_positions, layers=3)

        self.conv1_x = self._make_layer(block, 64, 1, 1)
        self.conv2 = BasicConv2D(input_channel=64,output_channel=128, strides=2)
        self.conv2_x = self._make_layer(block, 128, 1, 1)
        self.conv3 = BasicConv2D(input_channel=128,output_channel=256, strides=2)
        self.conv3_x = self._make_layer(block, 256, 1, 1)

        self.dropout = nn.Dropout(p=0.5)
        # self.conv3_x = self._make_layer(block, 64, num_block[1], 2)
        # self.conv4_x = self._make_layer(block, 256, num_block[2], 2)
        # self.conv5_x = self._make_layer(block, 512, num_block[3], 2)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(256 * depth_scale_size * width_scale_size, 16)
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
