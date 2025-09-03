%% 二值化函数
function inside = seg_easily(image,thread)
if nargin < 2
    thread = 0.7;  % 设置默认值
end
    % 转换为 uint8
    image = image - min(min(image));
    image_uint8 = uint8(double(image)*255 / max(image(:)));
%     figure()
%     imshow(image_uint8);
%     
    % 直方图均衡化
    equlized_image = histeq(image_uint8);
    
    % 阈值处理生成二值掩膜
%     level = min(1.3*graythresh(equlized_image),thread);
%     disp(level);

% 使用计算得到的阈值进行二值化
    inside = imbinarize(equlized_image, thread);
    inside = find_max_area_region(inside);
    inside =  imfill(inside, 'holes');
end