%% Prepare Data for SVM
% 处理64X128的图像，产生正样本参数和负样本参数，提供给SVM进行学习
% 一、图像以64X128为基准进行图像运算获得特征向量值
% 二、图像的WinSize(64, 128), BlockSize(16,16), CellSize(8,8), BlockStride(8,8),
% nbins(9)
% 三、图像特征向量值运算公式：
%             每个Cell有9个方向值
%             每个block有个2x2个Cell
%             每个Window有
% ((WinSize.x-BlockSize.x)/BlockStride.x+1)*((WinSize.y-BlockSize.y)/BlockStride.y+1)
%             既7x15个Block
%       所以64x128的图像应该有3780个特征向量

%% 清理系统参数，防止冲突
clear all;
close all;
%% 基本参数
cWinx = 64;cWiny = 128;
cBlkx = 16;cBlky = 16;
cCelx = 8;cCely = 8;
cSrdx = 8;cSrdy = 8;
%% 读取图像文件
%单个文件读取
%RdImage = imread('image/a0.jpg');
RdImage = imread('F:\opencv\sample\positive\a12.jpg');
%多个文件读取
%% 图像高斯滤波
%省略
%% 计算梯度
[Size,Ori] = ham_sobel2(RdImage);
%% 统计方向直方图
%descriptors = zeros((round((cWinx-cBlkx)/cSrdx)+1)*(round((cWiny-cBlky)/cSrdy)+1)*(cBlkx/cCellx)*(cBlky/cCelly),9);
%descriptors = ham_hist2(Size,Ori, cWinx, cWiny, cBlkx, cBlky, cCelx, cCely, cSrdx, cSrdy);
ham_hist3(Size, Ori,zeros(64,128), zeros(16,16), zeros(8,8), zeros(8,8), 9);
%% Normalize直方图
            
















