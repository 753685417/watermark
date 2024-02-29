%% Harris 特征点检测 具有随机性 
%针对旋转后的影像，选取特征点，建立向量，通过向量，计算旋转的角度
clc;
%---------------------------------------step3--------------------------------------------

emb_img_rotate1 = watermarked_img_rotate;

%提取rgb三维，对三维图像进行角点检测

%旋转后图像
R1=emb_img_rotate1(:,:,1);G1=emb_img_rotate1(:,:,2);B1=emb_img_rotate1(:,:,3);
emb_img_rotate_gray = 0.3*R1+0.59*G1+0.11*B1; 

%原始图像 
img_gray = 0.3*original_img_crop(:,:,1)+0.59*original_img_crop(:,:,2)+0.11*original_img_crop(:,:,3);
%% 原始图像检测
psf=fspecial('gaussian',[5,5],1); %高斯低通滤波
Ix=filter2([-1,0,1],img_gray);
Iy=filter2([-1,0,1]',img_gray);
Ix2=filter2(psf,Ix.^2);
Iy2=filter2(psf,Iy.^2);
Ixy=filter2(psf,Ix.*Iy);

[m,n]=size(img_gray);
R=zeros(m,n);%创建一个与灰度图相同的全0数组
max=0;
for i=1:m
    for j=1:n
        M=[Ix2(i,j),Ixy(i,j); Ixy(i,j),Iy2(i,j)];
        R(i,j)=det(M)-0.04*(trace(M))^2;
        if R(i,j)>max
            max=R(i,j);
        end
    end
end

thresh=0.95; %阈值可调
tmp=zeros(m,n);
neighbours=[-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
for i=2:m-1
    for j=2:n-1
        if R(i,j)>thresh*max
            for k=1:8
                if R(i,j)<R(i+neighbours(k,1),j+neighbours(k,2))
                    break;
                end
            end
            if k==8
                tmp(i,j)=1;
            end
        end
    end
end
figure(3);
subplot(121),imshow(img_gray),title('角点检测');
hold on; %保留原图像
for i=2:m-1
    for j=2:n-1
        if tmp(i,j)==1
            plot(j,i,'r.','MarkerSize', 40)
            
        end
    end
end
hold off;
%% 旋转后图像检测
psf=fspecial('gaussian',[5,5],1); %高斯低通滤波
Ix=filter2([-1,0,1],emb_img_rotate_gray);
Iy=filter2([-1,0,1]',emb_img_rotate_gray);
Ix2=filter2(psf,Ix.^2);
Iy2=filter2(psf,Iy.^2);
Ixy=filter2(psf,Ix.*Iy);

[m1,n1]=size(emb_img_rotate_gray);
R=zeros(m1,n1);%创建一个与灰度图相同的全0数组
max=0;
for i=1:m1
    for j=1:n1
        M=[Ix2(i,j),Ixy(i,j); Ixy(i,j),Iy2(i,j)];
        R(i,j)=det(M)-0.04*(trace(M))^2;
        if R(i,j)>max
            max=R(i,j);
        end
    end
end

tmp1=zeros(m1,n1);
neighbours=[-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
for i=2:m1-1
    for j=2:n1-1
        if R(i,j)>thresh*max
            for k=1:8
                if R(i,j)<R(i+neighbours(k,1),j+neighbours(k,2))
                    break;
                end
            end
            if k==8
                tmp1(i,j)=1;
            end
        end
    end
end

figure(3);
subplot(122),imshow(emb_img_rotate_gray),title('角点检测');
hold on; %保留原图像
for i=2:m1-1
    for j=2:n1-1
        if tmp1(i,j)==1
            plot(j,i,'g.','MarkerSize', 40)
            
        end
    end
end
hold off;
%% 特征点匹配
points_number=sum(tmp(:));points_number1=sum(tmp1(:));%特征点的数量

distence_total=zeros(points_number,points_number1);%特征点距离
position1=zeros(points_number,2);%匹配点位置
position2=zeros(points_number1,2);%匹配点位置
current_point=1;
for i=2:m-1 %原始图像的MN
    for j=2:n-1
        if tmp(i,j)==1           
             characteristic_region = img_gray(i-1:i+1,j-1:j+1);%原图像特征区域
             A = 1/9*(sum(characteristic_region(:)));%原图像周围像素均值
             A_ca=(characteristic_region-A).^2;
             S = 1/9*sum(A_ca(:));%原图像均方差     
             position1(current_point,1) = i;position1(current_point,2) = j;%
             
             current_point1=1;
             for k=2:m1-1
                 for l=2:n1-1
                     if tmp1(k,l)==1
                         characteristic_region1 = emb_img_rotate_gray(k-1:k+1,l-1:l+1);%攻击后图像特征区域
                         A1 = 1/9*(sum(characteristic_region1(:)));%攻击像周围像素均值
                         A_ca1=(characteristic_region1-A1).^2;
                         S1 = 1/9*sum(A_ca1(:));%原图像均方差
                         distence = 0.96*abs(A-A1)+1.17*abs(S-S1);%距离                        
                         distence_total(current_point,current_point1)=distence;
                         position2(current_point1,1) = k;position2(current_point1,2) = l;%
                         current_point1=current_point1+1;
                     end
                 end
             end
             current_point=1+current_point;
        end
    end
end
%
match=zeros(1,points_number);
for i=1:points_number
    Measurement = distence_total(i,:);
    [c,d] =min(Measurement);%c value min  d position
    match(1,i)=i;
end

%% 旋转
% 找到中心点，根据中心点进行
% angleDegrees1 =0;
% for i=1:points_number
%     point1 = position1(i,:);
%     point2 = position2(match(1,i),:);
%     vector1 = point1 -[m/2,n/2];          vector1_length = norm(vector1);
%     vector2 = point2 -[m1/2,n1/2];        vector2_length = norm(vector2);
%     %计算两个向量的点
    dotProduct = dot(vector1, vector2);
%     resize_percent = roundn(vector2_length/vector1_length,-2);
%     %模长
    normVec1 = norm(vector1);
    normVec2 = norm(vector2);
    angle = acos(dotProduct/(normVec1*normVec2));%弧度制
    angleDegrees = rad2deg(angle);%转角度
%     angleDegrees=roundn(angleDegrees,0);%保留
%     angleDegrees1 =angleDegrees1+angleDegrees;
% end
% rotateAngle=angleDegrees1/points_number;%ANGLE
% 
% %% 转回来
% emb_img_rotate_reversal = imrotate(emb_img_rotate1,-rotateAngle,'nearest','loose');%将图像顺时针旋转 * 度
% % imwrite(emb_img_rotate_reversal,"embedding-result/反转后图像.tif")
% 
% %裁剪由于旋转造成的多余部分
% [x_crop,y_crop] = find(emb_img_rotate_reversal > 0, 1, 'first');
% emb_img_rotate_reversal_crop = emb_img_rotate_reversal(x_crop:x_crop+m-1,y_crop:y_crop+n-1,:);
% imwrite(emb_img_rotate_reversal_crop,"Results/aa.tif")
% disp(strcat('旋转、','裁剪成功,被旋转角度为：',num2str(angleDegrees)));

%% 平移
% for i=1:points_number
%     point1 = position1(i,:);
%     point2 = position2(match(1,i),:);
%     circshift_x = point2(1,1)-point1(1,1);
%     circshift_y = point2(1,2)-point1(1,2);
% end
% [x_crop,y_crop] = find(J > 0, 1, 'first');
% emb_img_translate = J(x_crop:x_crop+m-1,y_crop:y_crop+n-1,:);

% imwrite(emb_img_translate,"embedding-result/emb_img_imtranslate1.tif")
% emb_img_circshift = circshift(emb_img,[1000,0]); 
%% 缩放
% if points_number>points_number1
%    [x,y]=find(distence_total==min(distence_total));
% end
% point1 = position1(x,:); vector1 = point1 -[m/2,n/2]; 
% point2 = position1(y,:); vector2 = point2 -[m1/2,n1/2];  
% resize_percent = roundn(vector2_length/vector1_length,-2);
% img_resize = imresize(emb_img_resize,1/resize_percent);



