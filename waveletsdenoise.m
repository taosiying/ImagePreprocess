%%基于小波变换的去噪函数：将图像分解为四个频带，分别为低频AA自带，HL子带，LH自带，HH子带，在AA子带保持不变，在其它子带分别用sobel算子计算边缘，在边缘区保持不变，其它区域通过软阈值变换，去除噪声
% 输入：noisy_img:带噪声图像，window_size1:对低频自带平滑窗口(已禁用，保持低频自带不变)，window_size2:对整体图像平滑窗口大小
% 使用示例：imgDen = waveletsdenoise(imgFlat，5，7);
%%
function denoised_img = waveletsdenoise(noisy_img,window_size1,window_size2)
if nargin < 2
        window_size1 = 5;  % 设置默认值n 
        window_size2 =7;
end 

noisy_img =  double(noisy_img);
diffused_img = noisy_img;
% diffused_img = anisotropic_diffusion(noisy_img, 20, 100, 0.25 ; %扩散去噪

% 定义 Sobel 滤波器
sobel_h = [1  2  1;   % 水平边缘检测
           0  0  0;
          -1 -2 -1];

sobel_v = [1  0 -1;   % 垂直边缘检测
           2  0 -2;
           1  0 -1];

sobel_d = [2  1  0;   % 对角边缘检测
           1  0 -1;
           0 -1 -2];
%% 子波包分解
% 小波包分解
wpt = wpdec2(diffused_img, 3, 'coif2');

% 不使用 wpjoin，直接使用原始树或最佳子树
best_tree = wpt; % 或者 best_tree = besttree(wpt);

% 获取叶节点
nodes = leaves(best_tree);
nodes = nodes(:);
if isempty(nodes) || any(nodes == 0)
    error('叶节点列表无效！nodes：%s', mat2str(nodes));
end

% 边缘检测与阈值处理
for n = nodes'
    coef = wpcoef(best_tree, n);    
    [depth, index] = getNodeInfo(n, best_tree); % 传递 best_tree 如果需要

    if mod(index, 4) == 0  % AA子带
        coef_new = coef;
        % coef_new = low_pass(coef, a); % 实现 low_pass 函数
%         coef_new = savitzky(coef,window_size1,3); % 实现 savitzky_1d 函数
%         coef_new = imguidedfilter(coef);   % 双边滤波

    elseif mod(index, 4) == 1  % HL
        edge_map = imfilter(coef, sobel_h, 'replicate') > 0;
        Thr = compute_threshold(coef); % 确保 compute_threshold 存在
        coef_new = coef;
        coef_new(~edge_map) = coef(~edge_map) .* (abs(coef(~edge_map)) > Thr);

    elseif mod(index, 4) == 2  % LH
        edge_map = imfilter(coef, sobel_v, 'replicate') > 0;
        Thr = compute_threshold(coef);
        coef_new = coef;
        coef_new(~edge_map) = coef(~edge_map) .* (abs(coef(~edge_map)) > Thr);

    elseif mod(index, 4) == 3  % HH
        edge_map = imfilter(coef, sobel_d, 'replicate') > 0;
        Thr = compute_threshold(coef);
        coef_new = coef;
        coef_new(~edge_map) = coef(~edge_map) .* (abs(coef(~edge_map)) > Thr);
    end

    best_tree = write(best_tree, 'data', n, coef_new);
end

% 重构图像
denoised_img = wprec2(best_tree);
denoised_img = savitzky(denoised_img,window_size2,1);
% figure()
% imagesc(denoised_img);
% figure()
% imagesc(noisy_img);

% 显示去噪结果
% figure;
% imshow(uint16(denoised_img)); title('Denoised Image');
%% 6. 性能评估
% 计算SNR
img = noisy_img;
SNR = 10 * log10(sum(img(:).^2) / sum((img(:) - denoised_img(:)).^2));

% 计算RMSE
RMSE = sqrt(mean((img(:) - denoised_img(:)).^2));

% 计算SSIM
ssim_val = ssim(uint8(denoised_img), uint8(img));

% 输出结果
fprintf('SNR: %.4f dB\n', SNR);
fprintf('RMSE: %.4f\n', RMSE);
fprintf('SSIM: %.4f\n', ssim_val);
% denoised_img = uint16(denoised_img);
end
function img_filtered = savitzky(img, window_size,poly_order)
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
end
function coef_new =  low_pass(coef,a)
 F = fft2(coef);
[h, w] = size(coef);
[X, Y] = meshgrid(1:w, 1:h);
cx = floor(w/2);
cy = floor(h/2);
radius = min(h, w) / a;
lp_filter = exp(-((X - cx).^2 + (Y - cy).^2) / (2 * (radius^2)));
F_filtered = F .* fftshift(lp_filter);
coef_new = real(ifft2(F_filtered));
end
function diffused_img = anisotropic_diffusion(img, iterations, K, lambda)
    diffused_img = img;
    [rows, cols] = size(img);
    
    for iter = 1:iterations
        % 计算梯度（北、南、东、西方向）
        grad_N = [diff(diffused_img, 1, 1); zeros(1, cols)];
        grad_S = -[zeros(1, cols); diff(diffused_img, 1, 1)];
        grad_E = [diff(diffused_img, 1, 2) zeros(rows, 1)];
        grad_W = -[zeros(rows, 1) diff(diffused_img, 1, 2)];
        
        % 计算导率系数
        c_N = exp(-(abs(grad_N)/K).^2);
        c_S = exp(-(abs(grad_S)/K).^2);
        c_E = exp(-(abs(grad_E)/K).^2);
        c_W = exp(-(abs(grad_W)/K).^2);
        
        % 更新图像
        diffused_img = diffused_img + lambda * (c_N .* grad_N + c_S .* grad_S + ...
            c_E .* grad_E + c_W .* grad_W);
    end
end
% 阈值计算函数
function Thr = compute_threshold(coef)
    [rows, cols] = size(coef);
    T = zeros(rows, cols);
    
    for i = 2:rows-1
        for j = 2:cols-1
            window = coef(i-1:i+1, j-1:j+1);
            cf = window(2, 2);
            ncf = [window(1,:), window(2,1), window(2,3), window(3,:)];
            
            m = mean(ncf);
            v = sqrt(mean((ncf - m).^2));
            d = cf - ncf;
            d_ave = mean(abs(d));
            dev = mean(abs(d - d_ave));
            T(i, j) = (d_ave + dev) * (m / (m + v));
        end
    end
    Thr = max(T(:));
end
% 辅助函数：获取节点深度和索引
function [depth, index] = getNodeInfo(node, best_tree)
    if ~isscalar(node) || node < 0
        error('节点索引必须是标量且非负！当前节点：%s', mat2str(node));
    end
    if node == 0
        depth = 0;
    elseif node <= 4
        depth = 1;
    elseif node <= 20
        depth = 2;
    else
        depth = 3; % 节点 21-84 为深度 3
    end
    index = mod(node - 21, 4); % 子带索引：0=AA, 1=HL, 2=LH, 3=HH
end