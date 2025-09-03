% 生成随机高斯曲面本经的函数
% 输入:img_size:图像尺寸，A_background：背景圆幅度，A_gaussian:高斯曲面幅度，num_gaussian:高斯波包个数
% 输出：与img_size相同的图像矩阵
% 使用示例：
% A_background = max(imgRaw(:)) * 0.3;
% A_gaussian = A_background;
% num_gaussian = 3;
% img_size = size(imgRaw,1);
% background = generate_background(img_size,A_background,A_gaussian,num_gaussian);
function result = generate_background(img_size,A_background,A_gaussian,num_gaussian)
% 设置参数
% img_size = 500;          % 图像尺寸
% 
% A_background = 200;
% A_gaussian = 50;
% num_gaussian = 3;
% 生成坐标网格
circle_radius = 0.6 * img_size;
[x, y] = meshgrid(1:img_size, 1:img_size);
center = [img_size/2, img_size/2]; % 中心位置
% 1. 生成大圆形背景（使用三角函数创建纹理）
circle_mask = ((x - center(1)).^2 + (y - center(2)).^2) <= circle_radius^2;
circle_mask = double(circle_mask) .* A_background;

background = zeros([img_size,img_size]) + circle_mask;

% 2. 生成随机高斯曲面
gaussian_surface = A_gaussian * generate_multiple_gaussians(x, y, img_size, num_gaussian);
gaussian_surface(background == 0) = 0;

% 3. 叠加背景和高斯曲面
result = background + gaussian_surface;

% 5. 3D可视化
figure;
imagesc(result);colorbar;

end


function gaussian_surface = generate_multiple_gaussians(x, y, N, num_gaussians)
    % 生成多个随机高斯曲面
    % 输入:
    %   x, y - 坐标网格
    %   N - 图像尺寸
    %   num_gaussians - 高斯曲面数量
    % 输出:
    %   gaussian_surface - 叠加后的高斯曲面
    
    gaussian_surface = zeros(size(x));
    center = [N/2, N/2]; % 图像中心
    
    % 高斯曲面参数范围
    amplitude_range = [1, 2.0];    % 振幅范围
    sigma_range = [N*0.3, N*0.5];         % 标准差范围
    center_spread = 0.9;             % 中心分布范围比例
    
    for i = 1:num_gaussians
        % 随机生成高斯曲面参数
        A = amplitude_range(1) + diff(amplitude_range) * rand();
        
        % 随机中心位置（70%概率在中心区域，30%概率在边缘区域）
        if rand() < 0.7
            % 在中心区域随机分布
            angle = 2 * pi * rand();
            max_distance = (N/2) * center_spread * (0.3 + 0.7 * rand());
            x0 = center(1) + max_distance * cos(angle);
            y0 = center(2) + max_distance * sin(angle);
        else
            % 在图像边缘区域随机分布
            edge_margin = N * 0.1; % 边缘留10%的边界
            if rand() < 0.5
                % 在左右边缘
                x0 = edge_margin + (N - 2*edge_margin) * rand();
                y0 = edge_margin;
            else
                % 在上下边缘
                x0 = N - edge_margin;
                y0 = edge_margin + (N - 2*edge_margin) * rand();
            end
        end
        
        % 确保中心在图像范围内
        x0 = max(1, min(N, x0));
        y0 = max(1, min(N, y0));
        
        % 随机方差（可以各向异性）
        sigma_x = sigma_range(1) + diff(sigma_range) * rand();
        sigma_y = sigma_range(1) + diff(sigma_range) * rand();
        
        % 50%的概率让高斯曲面更圆，50%的概率更椭圆
        if rand() < 0.5
            sigma_y = sigma_x * (0.8 + 0.4 * rand()); % 稍微椭圆
        else
            % 保持各向异性
        end
        
        % 生成高斯曲面
        gaussian = A * exp(-((x - x0).^2/(2*sigma_x^2) + (y - y0).^2/(2*sigma_y^2)));
        gaussian_surface = gaussian_surface + gaussian;
        
        fprintf('高斯曲面 %d: A=%.2f, 中心=(%.1f,%.1f), σx=%.1f, σy=%.1f\n', ...
                i, A, x0, y0, sigma_x, sigma_y);
    end
    
    % 可选：对高斯曲面进行归一化或调整整体强度
    max_val = max(gaussian_surface(:));
    if max_val > 0
        gaussian_surface = gaussian_surface / max_val * (amplitude_range(1) + diff(amplitude_range));
    end
end