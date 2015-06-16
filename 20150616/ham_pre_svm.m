%% Prepare Data for SVM
% ����64X128��ͼ�񣬲��������������͸������������ṩ��SVM����ѧϰ
% һ��ͼ����64X128Ϊ��׼����ͼ����������������ֵ
% ����ͼ���WinSize(64, 128), BlockSize(16,16), CellSize(8,8), BlockStride(8,8),
% nbins(9)
% ����ͼ����������ֵ���㹫ʽ��
%             ÿ��Cell��9������ֵ
%             ÿ��block�и�2x2��Cell
%             ÿ��Window��
% ((WinSize.x-BlockSize.x)/BlockStride.x+1)*((WinSize.y-BlockSize.y)/BlockStride.y+1)
%             ��7x15��Block
%       ����64x128��ͼ��Ӧ����3780����������

%% ����ϵͳ��������ֹ��ͻ
clear all;
close all;
%% ��������
cWinx = 64;cWiny = 128;
cBlkx = 16;cBlky = 16;
cCelx = 8;cCely = 8;
cSrdx = 8;cSrdy = 8;
%% ��ȡͼ���ļ�
%�����ļ���ȡ
%RdImage = imread('image/a0.jpg');
RdImage = imread('F:\opencv\sample\positive\a12.jpg');
%����ļ���ȡ
%% ͼ���˹�˲�
%ʡ��
%% �����ݶ�
[Size,Ori] = ham_sobel2(RdImage);
%% ͳ�Ʒ���ֱ��ͼ
%descriptors = zeros((round((cWinx-cBlkx)/cSrdx)+1)*(round((cWiny-cBlky)/cSrdy)+1)*(cBlkx/cCellx)*(cBlky/cCelly),9);
%descriptors = ham_hist2(Size,Ori, cWinx, cWiny, cBlkx, cBlky, cCelx, cCely, cSrdx, cSrdy);
ham_hist3(Size, Ori,zeros(64,128), zeros(16,16), zeros(8,8), zeros(8,8), 9);
%% Normalizeֱ��ͼ
            















