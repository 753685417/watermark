clc;
%% 攻击
%-------------step2----------------------------------

% 噪声攻击
% watermarked_img_noise=imnoise(watermarked_img,'salt & pepper',0.001);%添加方差为0.1的乘性噪声 speckle salt & pepper
% watermarked_img_noise=imnoise(watermarked_img,'gaussian',0,0.002);
watermarked_img_noise=imnoise(watermarked_img,'speckle',0.002);

%旋转
% watermarked_img_rotate = imrotate(watermarked_img,40,'nearest','loose');%将图像顺时针旋转*度

qu=floor((4096*4096*0.1)^(1/2));
%裁剪
img_crop = original_img;
% img_crop(4096-qu+1:4096,1:qu,:)=0;%左下
% img_crop(1:qu,4096-qu+1:4096,:)=0;%右上
% img_crop(4096-qu+1:4096,4096-qu+1:4096,:)=0;%右下
% img_crop(1000:1000+qu,1000:1000+qu,:)=0;%中间
img_crop(1:qu,1:qu,:)=0;%左上

%缩放
watermarked_img_resize = imresize(watermarked_img,0.9); %0.46
img_resize = imresize(watermarked_img_resize,[4096 4096]);

% %平移2
% watermarked_img_circshift = circshift(watermarked_img,[1500,1000]);
%平移
% J = imtranslate(watermarked_img,[0,100],'FillValues',0,'OutputView','full');

%低通滤波
% Iblur = imgaussfilt(watermarked_img,2.5); %标准差   苏丹2.4  %1 1.5 %subset1 1.9 %sh 3.0

%均值滤波
% m = 8;
% h = fspecial('average',[m m]);
% averageImageAttacked = imfilter(watermarked_img,h,'replicate');

%中值滤波
% m=6;
% % medianattack_img
% medianattack_img(:,:,1)=medfilt2(watermarked_img(:,:,1),[m,m]);
% medianattack_img(:,:,2)=medfilt2(watermarked_img(:,:,2),[m,m]);
% medianattack_img(:,:,3)=medfilt2(watermarked_img(:,:,3),[m,m]);
% medianattack_img(:,:,4)=medfilt2(watermarked_img(:,:,4),[m,m]);
% % medianfiltering = medfilt2(watermarked_img,[m,m]);


% imwrite(watermarked_img_noise,"embedding-result/noise.tif")%噪声
% imwrite(watermarked_img_rotate,"Results/旋转后.tif")%旋转
% imwrite(img_crop,"Results/watermarked_img_crop.tif")%裁剪

R=[0,-1;1,0;0,6000];
geotiffwrite('Results\img.tif', img_crop, R, 'CoordRefSysCode', 21417);

% imwrite(J,"embedding-result/watermarked_img_imtranslate.tif")%平移
%  imwrite(watermarked_img_circshift,"embedding-result/watermarked_img_circshift.tif")%平移2
 
%a = PSNR(watermarked_img_noise,original_img_crop);
% imwrite(averageImageAttacked,"embedding-result/averageImageAttacked.tif")%均值滤波
disp('攻击结束');