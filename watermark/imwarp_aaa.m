% clc,clear,close all;
% 
% watermarking_img = imread('QRdata\苏州科技大学.png');%读取遥感影像和水印文件
% watermarking_img = im2gray(watermarking_img);%水印转为灰度图，同时降维

A = [1 0.002 0;   
     0 1 0; 
     0 0 1];
tform = affinetform2d(A);

J = imwarp(watermarked_img,tform);

k = imresize(J,[4096,4096]);
