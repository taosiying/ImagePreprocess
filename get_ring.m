%% 计算半影区主函数
function [ring_mask,outside_mask,inside_mask] = get_ring(image,dr, xc, yc, a, b, theta)
    
    % 计算拟合参数
    num_sample = 500;
    [~, a_b_inside, a_b_outside] = cal_fit(image, xc, yc, a, b, theta, num_sample);
    
    % 获取图像尺寸
    [rows, cols] = size(image);
    
    % 创建内外椭圆掩膜
    [X, Y] = meshgrid(1:cols, 1:rows);
    x_rot = cos(theta) * (X - xc) + sin(theta) * (Y - yc);
    y_rot = -sin(theta) * (X - xc) + cos(theta) * (Y - yc);
    
    % 外部椭圆掩膜
    outside_mask = uint8((x_rot / (a_b_outside(1) + dr)).^2 + (y_rot / (a_b_outside(2) + dr)).^2 <= 1) * 255;
    
    % 内部椭圆掩膜
    inside_mask = uint8((x_rot / (a_b_inside(1) - dr)).^2 + (y_rot / (a_b_inside(2) - dr)).^2 <= 1) * 255;
    
    % 生成环形掩膜
    ring_mask = bitand(outside_mask, bitcmp(inside_mask, 'uint8'));

end

%% 线性拟合函数
  function [r_start, r_end] = linfit(data,show_result )
    % 分段线性拟合，检测显著下降区域
    x = 1:length(data);
    n_segments = 3;
    
    % 检查数据中的 NaN 或 Inf
    valid_idx = isfinite(data); % 排除 NaN 和 Inf
    if sum(valid_idx) < 2 % 至少需要2个有效点进行拟合
        warning('输入数据包含过多 NaN/Inf 或为空，返回默认值。');
        r_start = 1;
        r_end = length(data);
        return;
    end
    
    % 移除 NaN/Inf 值用于拟合
    x_clean = x(valid_idx);
    data_clean = data(valid_idx);
    
    % 使用 fitnlm 进行线性拟合（模拟 Python 的 pwlf）
    try
        mdl = fitnlm(x_clean', data_clean', 'y ~ b0 + b1*x', [mean(data_clean), 0]);
    catch e
        warning('fitnlm 失败：%s，返回默认值。', e.message);
        r_start = 1;
        r_end = length(data);
        return;
    end
    
    % 计算分段斜率
    slopes = [];
    breakpoints = round(linspace(1, length(data), n_segments + 1));
    
    for i = 1:(length(breakpoints) - 1)
        start = breakpoints(i);
        end_idx = breakpoints(i + 1);
        if end_idx > length(data) || start > length(data)
            continue; % 跳过无效分段
        end
        slope = (data(end_idx) - data(start)) / (end_idx - start);
        amplitude = data(start) - data(end_idx);
        slopes = [slopes; [start, end_idx, slope, amplitude]];
    end
    
    % 检查 slopes 是否为空
    if isempty(slopes)
        warning('无有效分段用于斜率计算，返回默认值。');
        r_start = 1;
        r_end = length(data);
        return;
    end
    
    % 找到绝对斜率最大的分段
    [~, max_index] = max(abs(slopes(:, 3)));
    r_start = max(1, slopes(max_index, 1) - 20);
    r_end = min(length(data), slopes(max_index, 2) + 20);
    
    % 如果需要，绘制结果
    if 1
        figure;
        plot(x, data, 'b-', 'DisplayName', '数据');
        hold on;
        xline(r_start, 'r--', 'DisplayName', '起点');
        xline(r_end, 'g--', 'DisplayName', '终点');
        xlabel('索引');
        ylabel('值');
        title('分段线性拟合');
        legend;
        grid on;
        drawnow;
    end
  end

%% 扩张缩小椭圆统计各个半径下灰度均值,计算内外半径的函数
function [grays, a_b_inside, a_b_outside] = cal_fit(image, xc, yc, a, b, theta, num_sample)
    % Compute grayscale means along scaled ellipses and find inner/outer boundaries
    scales = linspace(0.3, 1.5, num_sample);
    grays = zeros(1, num_sample);
    
    for i = 1:num_sample
        a1 = a * scales(i);
        b1 = b * scales(i);
        grays(i) = compute_ellipse_boundary_gray_mean(image, xc, yc, a1, b1, theta);
    end
    
    [start_idx, end_idx] = piecewise_linear_fit(grays);
%     fprintf("start idx:%d,end idx:%d",start_idx,end_idx)
%     fprintf("a_inside:%f,b_inside%f",a * scales(start_idx), b * scales(start_idx))
    a_b_inside = [a * scales(start_idx), b * scales(start_idx)];
    a_b_outside = [a * scales(end_idx), b * scales(end_idx)];
end

%% 计算某个椭圆下的灰度均值
function gray_mean = compute_ellipse_boundary_gray_mean(gray_img, xc, yc, a, b, theta)
    % Compute mean grayscale value along an ellipse boundary
    mask = zeros(size(gray_img), 'uint8');
    [X, Y] = meshgrid(1:size(gray_img, 2), 1:size(gray_img, 1));
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    x_rot = cos_theta * (X - xc) + sin_theta * (Y - yc);
    y_rot = -sin_theta * (X - xc) + cos_theta * (Y - yc);
    ellipse = abs((x_rot / a).^2 + (y_rot / b).^2 - 1) < 0.01; % Thin boundary
    boundary_pixels = gray_img(ellipse);
    if isempty(boundary_pixels)
        gray_mean = NaN;
    else
        gray_mean = mean(boundary_pixels);
    end
end

%% 分段线性拟合
function [start_idx, end_idx] =  piecewise_linear_fit(data)
    N = length(data);
    min_seg_len = floor(N / 10);  % 控制最小段长度

    best_slope = -inf;
    start_idx = 1;
    end_idx = N;
    max_slope = 0;

    % 枚举两断点 i, j
    for i = min_seg_len : (N - 2 * min_seg_len)
        for j = (i + min_seg_len) : (N - min_seg_len)
            seg1 = 1:i;
            seg2 = (i+1):j;
            seg3 = (j+1):N;

            k1 = ols_slope(seg1, data(seg1));
            k2 = ols_slope(seg2, data(seg2));
            k3 = ols_slope(seg3, data(seg3));

            slopes = [k1, k2, k3];
            segments = [1, i; i+1, j; j+1, N];

            [slope_val, idx] = max(abs(slopes));
            if slope_val > best_slope
                best_slope = slope_val;
                start_idx = segments(idx, 1);
                end_idx = segments(idx, 2);
                max_slope = slopes(idx);
            end
        end
    end
end
function k = ols_slope(x, y)
    x_mean = mean(x);
    y_mean = mean(y);
    num = sum((x - x_mean) .* (y - y_mean));
    den = sum((x - x_mean).^2);
    k = num / den;
end