%% 对噪声或其他攻击直接提取
clc;
extract_img = original_img_crop ;

final_image = zeros(32);
nc_value = 0;
for q=1:3
    for w =1:3
        caculate_place = extract_img((q-1)*512+1:q*512,(w-1)*512+1:w*512,:);

        dwt3_result12 = dwt3(caculate_place, wname);%嵌入结果、攻击结果
        lowpass_filtering12 = dwt3_result12.dec{1,1,1};
        dwt3_result22 = dwt3(lowpass_filtering12, wname);
        lowpass_filtering22 = dwt3_result22.dec{1,1,1};
        %dct
        dct_caculated2 = blockproc(lowpass_filtering22,[8 8],dct);
        %提取嵌入区域
        Embeddingplace_cell2 = cell(8, 8);
        for l = 1:8
            for m = 1:8
                Embeddingplace_cell2{l,m} = dct_caculated2((l-1)*block_width+1:(l-1)*block_width+2, (m-1)*block_width+1:(m-1)*block_width+2);
            end
        end
        %转换为矩阵
        Embeddingplace_new_return = cell2mat(Embeddingplace_cell2);
        %嵌入区域svd
        [Uz2,Sz_emb2,Vz2] = svd(Embeddingplace_new_return);
        Sw2 = (Sz_emb2-Sz)/alpha;

        cA_2 = Uw*Sw2*Vw';
        extraction_Wmdata = idwt2(cA_2,cH,cV,cD,'haar');
        %反置乱
        extraction_Wmdata = im2uint8(extraction_Wmdata);%转回uint8
        rearnoldimg = rearnold(extraction_Wmdata,3,5,10); 
        %将水印值重置为0、255
        rearnoldimg(rearnoldimg<125)= 0;rearnoldimg(rearnoldimg>=125)= 255;

        dnc = nc(rearnoldimg,watermarking_img);%函数nc
       

        if nc_value<dnc
            nc_value = dnc;
            final_image = rearnoldimg;
        end
    end
end
imwrite(final_image,"Results/extraction_Wmdata.png")

figure(2);
subplot(1,2,1), imshow(watermarking_img),title('original watermaring');
subplot(1,2,2), imshow(final_image),title(strcat('extract','  NC=',num2str(nc_value)));%串联
disp(nc_value);