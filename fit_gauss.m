%%
% 对图像image的mask区域的点按照sample_ratio比例采样，并拟合为num_gaussians个高斯函数组成的曲面
% image：输入图像数据
% mask：被拟合区域二值图像
% num_gaussians:高斯函数个数
% sample_ratio:采样频率
% max_iterations:最大迭代次数
% max_sigma：最大方差比例，限制方差为max_sigma*size(image,1)
% constrain_flat:是否限制mask以外区域尽量平坦
% fitted_surface:拟合曲面结果
%%

function [fitted_surface, gaussian_models, all_params, rmse] = fit_gauss(image, mask, num_gaussians, sample_ratio, max_iterations, max_sigma, constrain_flat)
% 添加 constrain_flat 参数来控制是否约束mask外区域平坦
if nargin < 7
    constrain_flat = false; % 默认不约束mask外区域平坦
end
if nargin < 6
    max_sigma = 2;
end

[rows, cols] = find(mask > 0);
intensities = double(image(sub2ind(size(image), rows, cols)));

num_points = numel(rows);
image_size = size(image,1);
if num_points == 0
    error('mask区域没有找到任何像素点');
end

% 随机采样
sample_count = round(num_points * sample_ratio);
sample_count = min(sample_count, num_points);
rand_indices = randperm(num_points, sample_count);

% 准备拟合数据
X = cols(rand_indices);
Y = rows(rand_indices);
Z = intensities(rand_indices);

% 初始化所有高斯函数的参数 - 使用mask内的数据
all_params = initializeParameters(X, Y, Z, num_gaussians, mask,image_size);

% 合并所有参数为一个向量
initial_params = all_params(:);

% 设置参数边界 - 限制在mask区域内，并设置最小方差
[lb, ub] = setParameterBounds(X, Y, Z, num_gaussians, mask, max_sigma,image_size);

% 拟合选项
options = optimset('Display', 'iter', 'MaxIter', max_iterations, 'TolFun', 1e-6, 'TolX', 1e-6);

% 如果需要约束mask外区域平坦
if constrain_flat
    % 获取mask外的采样点用于平坦性约束
    [out_rows, out_cols] = find(mask == 0);
    if ~isempty(out_rows)
        % 随机采样mask外的点
        out_num_points = numel(out_rows);
        out_sample_count = min(round(out_num_points * sample_ratio * 0.3), 1000); % 限制采样数量
        out_rand_indices = randperm(out_num_points, out_sample_count);
        
        X_out = out_cols(out_rand_indices);
        Y_out = out_rows(out_rand_indices);
        
        % 合并所有数据点
        X_all = [X; X_out];
        Y_all = [Y; Y_out];
        Z_all = [Z; zeros(size(X_out))]; % mask外期望值为平坦
        
        % 使用自定义目标函数，包含平坦性约束
        fitted_params = lsqcurvefit(@(p,xy) flatnessConstrainedMultiGaussian2D(p, xy, num_gaussians, length(X)), ...
                                   initial_params, [X_all, Y_all], Z_all, lb, ub, options);
    else
        % 如果没有mask外区域，使用普通拟合
        fitted_params = lsqcurvefit(@(p,xy) multiGaussian2D(p,xy,num_gaussians), ...
                                   initial_params, [X, Y], Z, lb, ub, options);
    end
else
    % 原始拟合方法
    fitted_params = lsqcurvefit(@(p,xy) multiGaussian2D(p,xy,num_gaussians), ...
                               initial_params, [X, Y], Z, lb, ub, options);
end

% 重新组织参数
all_params = reshape(fitted_params, num_gaussians, 7);

% 计算拟合质量 - 只在mask区域内计算
Z_fit = multiGaussian2D(fitted_params, [X, Y], num_gaussians);
rmse = sqrt(mean((Z - Z_fit).^2));

% 创建各个高斯模型函数
gaussian_models = cell(num_gaussians, 1);
for i = 1:num_gaussians
    params = all_params(i, :);
    gaussian_models{i} = @(x,y) singleGaussian2D(params, [x(:), y(:)]);
end

% 生成与原始图像相同尺度的拟合曲面
fitted_surface = generateFittedSurface(image, all_params, num_gaussians);

% 显示拟合结果
displayMultiGaussianResults(image, mask, X, Y, Z, fitted_surface, all_params, rmse);
end

% 修改：简化平坦性约束函数，正确处理输入数据
function Z = flatnessConstrainedMultiGaussian2D(params, XY, num_gaussians, n_in_points)
% 多个二维高斯函数的叠加，包含mask外平坦性约束
% XY 应该是一个 N×2 的矩阵，第一列是X坐标，第二列是Y坐标

% 检查输入数据格式
if size(XY, 2) ~= 2
    error('XY 应该是一个 N×2 的矩阵，包含X和Y坐标');
end

params = reshape(params, num_gaussians, 7);
Z_pred = zeros(size(XY,1), 1);

for i = 1:num_gaussians
    Z_pred = Z_pred + singleGaussian2D(params(i,:), XY);
end

% 分离mask内外的预测值
n_total = size(XY, 1);
n_out_points = n_total - n_in_points;

% 计算mask内的拟合值
Z_in_pred = Z_pred(1:n_in_points);

% 计算mask外的平坦性约束
if n_out_points > 0
    Z_out_pred = Z_pred(n_in_points+1:end);
    
    % 计算mask外区域的平坦性（使用标准差作为平坦性度量）
    out_flatness = std(Z_out_pred) * 0.1; % 平坦性误差，权重较低
    
    % 组合结果：mask内拟合值 + mask外平坦性约束
    Z = [Z_in_pred; out_flatness * ones(n_out_points, 1)];
else
    % 如果没有mask外点，只返回mask内预测值
    Z = Z_in_pred;
end
end


function Z = multiGaussian2D(params, XY, num_gaussians)
    % 多个二维高斯函数的叠加
    params = reshape(params, num_gaussians, 7);
    Z = zeros(size(XY,1), 1);
    
    for i = 1:num_gaussians
        Z = Z + singleGaussian2D(params(i,:), XY);
    end
end

function Z = singleGaussian2D(params, XY)
    % 单个二维旋转高斯函数
    A = params(1);
    x0 = params(2);
    y0 = params(3);
    sigma_x = params(4);
    sigma_y = params(5);
    theta = params(6);
    offset = params(7);
    
    x = XY(:,1);
    y = XY(:,2);
    
    % 坐标旋转
    a = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
    b = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
    c = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
    
    % 高斯函数
    exponent = a*(x-x0).^2 + 2*b*(x-x0).*(y-y0) + c*(y-y0).^2;
    Z = A * exp(-exponent) + offset;
end

function params = initializeParameters(X, Y, Z, num_gaussians, mask,image_size)
    % 初始化多个高斯函数的参数 - 只在mask区域内初始化，并设置合理的方差
    params = zeros(num_gaussians, 7);
    
    max_intensity = max(Z);
    min_intensity = min(Z);
    
    % 获取mask的边界框
    [mask_rows, mask_cols] = find(mask > 0);
    min_x = min(mask_cols);
    max_x = max(mask_cols);
    min_y = min(mask_rows);
    max_y = max(mask_rows);
    
    range_x = max_x - min_x;
    range_y = max_y - min_y;
    
    % 计算mask区域的面积和特征尺寸
%     mask_area = sum(mask(:) > 0);
%     characteristic_size = sqrt(mask_area / num_gaussians);
    
    for i = 1:num_gaussians
        % 在mask区域内随机分布中心点
        x0 = min_x + range_x * rand();
        y0 = min_y + range_y * rand();
        
        % 确保初始中心点在mask内
        attempts = 0;
        while ~isPointInMask(x0, y0, mask) && attempts < 100
            x0 = min_x + range_x * rand();
            y0 = min_y + range_y * rand();
            attempts = attempts + 1;
        end
        
        % 设置合理的初始方差，避免点状曲面
        % 方差与特征尺寸成正比，确保高斯分布有足够的宽度
        base_sigma = image_size * 0.3;  % 基础方差
        sigma_factor = 0.8 + 0.4 * rand();       % 随机因子增加多样性
        
        sigma_x_init = base_sigma * sigma_factor;
        sigma_y_init = base_sigma * sigma_factor * (0.8 + 0.4 * rand());
        
        params(i, :) = [
            (max_intensity - min_intensity) / num_gaussians * (0.8 + 0.4*rand()),  % 幅度
            x0,                        % x中心
            y0,                        % y中心
            sigma_x_init,              % x方向标准差
            sigma_y_init,              % y方向标准差
            pi * (rand() - 0.5),       % 旋转角度
            min_intensity * (0.2 + 0.6*rand())  % 背景偏移
        ];
    end
end

function [lb, ub] = setParameterBounds(X, Y, Z, num_gaussians, mask,max_sigma,image_size)
    % 设置参数边界 - 限制在mask区域内，并设置方差的最小最大值
    max_intensity = max(Z);
    min_intensity = min(Z);
    max_offset = max(abs(max_intensity),abs(min_intensity));
    
    % 获取mask的边界框
    [mask_rows, mask_cols] = find(mask > 0);
    min_x = min(mask_cols);
    max_x = max(mask_cols);
    min_y = min(mask_rows);
    max_y = max(mask_rows);
    
    range_x = max_x - min_x;
    range_y = max_y - min_y;
    
    % 计算合理的方差范围
%     mask_area = sum(mask(:) > 0);
%     characteristic_size = sqrt(mask_area);
    characteristic_size = image_size;
    % 设置方差的最小和最大值
    min_sigma = max(2, characteristic_size * 0.05);  % 最小方差，避免点状
    max_sigma = characteristic_size * max_sigma;           % 最大方差，避免过平滑
    
    % 单个高斯函数的边界
    single_lb = [min(0,min_intensity), min_x, min_y, min_sigma, min_sigma, -pi, -max_offset];
    single_ub = [(max_intensity - min_intensity), max_x, max_y, max_sigma, max_sigma, pi, max_offset];
    
    % 扩展到所有高斯函数
    lb = repmat(single_lb, num_gaussians, 1);
    ub = repmat(single_ub, num_gaussians, 1);
    
    lb = lb(:)';
    ub = ub(:)';
end

function in_mask = isPointInMask(x, y, mask)
    % 检查点是否在mask内
    x_round = round(x);
    y_round = round(y);
    
    if x_round < 1 || x_round > size(mask, 2) || y_round < 1 || y_round > size(mask, 1)
        in_mask = false;
    else
        in_mask = mask(y_round, x_round) > 0;
    end
end

function fitted_surface = generateFittedSurface(image, all_params, num_gaussians)
    % 生成与原始图像相同尺度的拟合曲面
    [m, n] = size(image);
    [X_grid, Y_grid] = meshgrid(1:n, 1:m);
    XY = [X_grid(:), Y_grid(:)];
    
    % 计算所有高斯函数的叠加
    Z_fit = zeros(size(XY,1), 1);
    for i = 1:num_gaussians
        Z_fit = Z_fit + singleGaussian2D(all_params(i,:), XY);
    end
    
    fitted_surface = reshape(Z_fit, m, n);
end

function displayMultiGaussianResults(image, mask, X, Y, Z, fitted_surface, all_params, rmse)
    % 显示多高斯拟合结果
    
    figure('Position', [100, 100, 1200, 800]);
    
    % 原始图像和mask
    subplot(2,3,1);
    imshow(image, []);
    hold on;
    contour(mask, 'r', 'LineWidth', 1);
    title('原始图像和mask区域');
    colorbar;
    
    % 采样点
    subplot(2,3,2);
    scatter3(X, Y, Z, 10, Z, 'filled');
    colorbar;
    title('采样点强度分布');
    xlabel('X'); ylabel('Y'); zlabel('强度');
    
    % 拟合曲面
    subplot(2,3,3);
    imagesc(fitted_surface);colorbar;
    title(sprintf('多高斯拟合曲面 (RMSE=%.3f)', rmse));
    colorbar;
    
    % 残差图
    subplot(2,3,4);
    [x_grid, y_grid] = meshgrid(linspace(min(X), max(X), 50), linspace(min(Y), max(Y), 50));
    xy_sample = [x_grid(:), y_grid(:)];
    z_fit_sample = zeros(size(xy_sample,1), 1);
    for i = 1:size(all_params,1)
        z_fit_sample = z_fit_sample + singleGaussian2D(all_params(i,:), xy_sample);
    end
    z_fit_sample = reshape(z_fit_sample, size(x_grid));
    
    surf(x_grid, y_grid, z_fit_sample, 'EdgeColor', 'none');
    hold on;
    scatter3(X, Y, Z, 20, 'r', 'filled');
    title('三维拟合曲面');
    xlabel('X'); ylabel('Y'); zlabel('强度');
    
    % 各个高斯分量
    subplot(2,3,5);
    [m, n] = size(image);
    individual_surfaces = zeros(m, n, size(all_params,1));
    for i = 1:size(all_params,1)
        individual_surfaces(:,:,i) = generateFittedSurface(image, all_params(i,:), 1);
    end
    
    montage(individual_surfaces, 'DisplayRange', []);
    title('各个高斯分量');
    
    % 参数显示
    subplot(2,3,6);
    axis off;
    text(0.1, 0.9, sprintf('多高斯拟合参数 (%d个高斯):', size(all_params,1)), 'FontSize', 9);
    
    for i = 1:size(all_params,1)
        params = all_params(i,:);
        y_pos = 0.8 - (i-1)*0.12;
        text(0.1, y_pos, sprintf('高斯%d: A=%.1f, (x0,y0)=(%.1f,%.1f)', i, params(1), params(2), params(3)), 'FontSize', 8);
        text(0.1, y_pos-0.04, sprintf('     σ=(%.1f,%.1f), θ=%.2f, offset=%.1f', params(4), params(5), params(6), params(7)), 'FontSize', 8);
    end
    
    text(0.1, 0.1, sprintf('总RMSE = %.3f', rmse), 'FontSize', 10, 'Color', 'red');
end