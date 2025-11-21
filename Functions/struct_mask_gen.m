function [mask] = struct_mask_gen(ori_img, mask_type)
    if strcmp(mask_type, 'number')
        % 在原图上写数字
        loc = int32([0 0;0 170;0 80]);
        Y = insertText(ori_img,loc,'1234567','FontSize',50, 'BoxOpacity',0);
        Y = rgb2gray(Y);%将图像灰度化
        Y(Y ~= 1) = 0;
        % 生成mask_number
        mask_number = ones(size(ori_img));
        mask_number(:,:,1) = Y;
        mask_number(:,:,2) = Y;
        mask_number(:,:,3) = Y;
        mask = mask_number;
        return
    elseif strcmp(mask_type,'strip')
        [n1,n2,~] = size(ori_img);
        X = ones( n1, n2 );
%         a = [2:5, 44:46, 74:77, 80:81, 102:103, 134:135, 165:166, 183:184, 223:226, 246:247];    % 2:5 从第二行开始，宽度为4的条带    % 随机产生
%         b = [16:17, 26:30, 44:45, 54:55, 86:87, 136:137, 139:141, 200:201, 238:241, 249:250];    % 26:30 从第26列开始，宽度为5的条带
%         X(a,:) = 0;
%         X(:,b) = 0;
% 块缺失
        X(10:30,80:100) = 0;
        X(110:130,20:50) = 0;
        X(110:130,200:220) = 0;
        X(200:220,120:140) = 0;
        % 生成mask_strip
        mask_strip = ones(size(ori_img));
        mask_strip(:,:,1) = X;
        mask_strip(:,:,2) = X;
        mask_strip(:,:,3) = X;
        mask = mask_strip;
        return
    end
end

