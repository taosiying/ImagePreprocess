
% SG滤波平滑函数
% img：输入图像
% window_size:平滑窗口大小

function img_filtered = savitzky_1d(img, window_size)
    % 默认窗口大小为 51，如果未提供
    if nargin < 2
        window_size = 51;
    end
    
    % 多项式阶数
    poly_order = 1;
    
    % 获取图像尺寸
    [rows, cols] = size(img);
    
    % 初始化输出图像
    img_filtered = zeros(rows, cols, 'single');
    
    % 对每一列应用 Savitzky-Golay 滤波
    for j = 1:cols
        img_filtered(:, j) = sgolayfilt(double(img(:, j)), poly_order, window_size);
    end
    
    % 对每一行应用 Savitzky-Golay 滤波
    for i = 1:rows
        img_filtered(i, :) = sgolayfilt(double(img_filtered(i, :)), poly_order, window_size);
    end
    
    % 限制输出值在 [0, 65535] 范围内并转换为 uint16
% %     img_filtered = max(min(img_filtered, 65535), 0);
%     img_filtered = uint16(img_filtered);
end
