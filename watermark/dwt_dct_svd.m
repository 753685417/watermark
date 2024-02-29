
%---------嵌入----------step1-------------------------
clc,clear,close all;
%读取
original_img = imread('E:\data\1\subset.tif');
original_img = original_img(:,:,1:4);
%nanjing (3001:4536,3001:4536,:)
watermarking_img = imread('E:\data\watermark\苏州科技大学.png');%读取遥感影像和水印文件
watermarking_img = im2gray(watermarking_img);%水印转为灰度图，同时降维
% watermarking_img =imresize(watermarking_img,[]);

%对水印进行置乱
arnoldImg = arnold(watermarking_img,3,5,10);%a 3,b 5,n 10
arnoldImg = im2double(arnoldImg);%将水印值归到0-1间

original_img_crop = original_img(3001:4536,3001:4536,:);%计算区域 %256的整数倍 1536 2001:3536,2001:3536
% original_img_crop = original_img(1:1536,1:1536,:);
original_img_crop = mat2gray(original_img_crop);

%将水印dwt分解
[cA,cH,cV,cD] = dwt2(arnoldImg,'haar');%可实现降低嵌入量
%将低频分量，进行svd分解
[Uw,Sw,Vw] = svd(cA);%对HH变奇异值分解 Sw待用

%离散余弦函数
T = dctmtx(8);
dct = @(block_struct) T * block_struct.data * T';
wname = 'haar';%设定haar算子
emb_img = original_img_crop;

for q=1:3
    for w =1:3
        original_img_cacualte = original_img_crop((q-1)*512+1:q*512,(w-1)*512+1:w*512,:);

        %对遥感影像进行8*8 dwt变换,并进行dct
       
        dwt3_result1 = dwt3(original_img_cacualte, wname);
        lowpass_filtering1 = dwt3_result1.dec{1,1,1};
        dwt3_result2 = dwt3(lowpass_filtering1, wname);
        lowpass_filtering2 = dwt3_result2.dec{1,1,1};

        dct_caculated = blockproc(lowpass_filtering2,[8 8],dct);

        [length,width] = size(dct_caculated);
        block_width = length/8;
        Embeddingplace_cell = cell(8, 8);
        for i = 1:8
            for j = 1:8
                Embeddingplace_cell{i,j} = dct_caculated((i-1)*block_width+1:(i-1)*block_width+2, (j-1)*block_width+1:(j-1)*block_width+2);
            end
        end

        Embeddingplace = cell2mat(Embeddingplace_cell);%转换为矩阵
        %嵌入区域svd
        [Uz,Sz,Vz] = svd(Embeddingplace);

        alpha = 0.1; %嵌入强度--------------------------------------------------------
        Sz_emb = Sz + alpha*Sw; %嵌入加法
        % Sz_emb = Sz*(1+alpha*Sw); %乘法
        Embeddingplace_dct_new = Uz*Sz_emb*Vz';

        %获得dct变换后的矩阵
        dct_caculated_return = dct_caculated;
        for i = 1:8
            for j = 1:8
                dct_caculated_return((i-1)*block_width+1:(i-1)*block_width+2, (j-1)*block_width+1:(j-1)*block_width+2) = Embeddingplace_dct_new((i-1)*2+1:i*2, (j-1)*2+1:j*2);
            end
        end
        %将恢复的数据反余弦
        invdct = @(block_struct) T' * block_struct.data * T;
        dct_caculated_invdct = blockproc(dct_caculated_return,[8 8],invdct);

        lowpass_filtering2_return = dct_caculated_invdct;
        dwt3_result2.dec{1,1,1} = lowpass_filtering2_return;
        lowpass_filtering1_return = idwt3(dwt3_result2);
        dwt3_result1.dec{1,1,1} = lowpass_filtering1_return;
        emb_img1 =  idwt3(dwt3_result1);

        emb_img((q-1)*512+1:q*512,(w-1)*512+1:w*512,:) = emb_img1;
    end
end




%PSNR
dPSNR = psnr(emb_img,original_img_crop);
%ssim
%ssim_value = SSIM(emb_img,original_img_cacualte);

disp(dPSNR);disp('嵌入完成');
% imwrite(emb_img,"embedding-result/emb_img.tif")


