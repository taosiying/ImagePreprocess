function  RadialProfile_cmp(imgRaw,imgNew,x0,y0,offset,iplot,string,imask)
    %% 主脚本：8 射线强度对比
% 读取两幅图像并统一为 double
if nargin < 8
    imask = false;
end
im_origin = im2double(imgRaw);          % 原始图
im_proc   = im2double(imgNew);          % 处理后图

% 提取 8 条射线强度
intens_orig = extract8Lines(im_origin, 2048,x0,y0,iplot);
intens_proc= extract8Lines(im_proc,   2048,x0,y0,iplot);
if imask
   binaryraw = seg_easily(imgRaw);
   [xc,yc,a,b,theta] = fit_ellipse_lsp(binaryraw);
   ring_mask_raw = get_ring(imgRaw,0,xc,yc,a,b,theta);
   binarynew = seg_easily(imgNew);
   [xc,yc,a,b,theta] = fit_ellipse_lsp(binarynew);
   ring_mask_new = get_ring(imgRaw,0,xc,yc,a,b,theta);
   intens_mask_raw = max(imgRaw(:)) * extract8Lines(ring_mask_raw,   2048,x0,y0,iplot);
   intens_mask_new = max(imgNew(:)) * extract8Lines(ring_mask_new,   2048,x0,y0,iplot);
end
% 绘制 8 幅子图
figure('Position', [100 100 1500 900]);
offset = offset;                            % 垂直平移量，便于区分
angles = (0:22.5:157.5);                   % 8 个角度（与函数一致）
% 获取 figure 位置
for k = 1:8
    subplot(2,4,k);
    plot(intens_orig(k,:), 'r', 'LineWidth', 1.3); hold on;
    plot(intens_proc(k,:) + offset, 'b', 'LineWidth', 1.3);
    if imask
        plot(intens_mask_raw(k,:), 'r', 'LineWidth', 1);hold on;
        plot(intens_mask_new(k,:), 'b', 'LineWidth', 1);hold on;
    end
    grid on;
    title(string);
%     title(sprintf([string,' :%.1f°'], angles(k)));
    legend('Original', 'Processed', 'Location', 'best');
end
end

function intens = extract8Lines(I, len,x0,y0,iplot,interpMethod)
% extract8Lines  以图像中心为原点，沿 8 条均匀角度射线提取强度
%
% 输入
%   I            : 单幅图像（灰度或 RGB，任意尺寸）
%   len          : 射线半长（像素）
%   interpMethod : 插值方法，默认 'bilinear'
%
% 输出
%   intens       : 8×(2*len+1) 矩阵，每行一条射线强度
%   figHandle    : 叠加 8 条线的图像句柄（可忽略）

    if nargin < 6, interpMethod = 'bilinear'; end
    I = im2double(I);
    [H, W, ~] = size(I);
    %xc = (W+1)/2;   yc = (H+1)/2;
    xc=x0;yc=y0;
    
    theta = (0:22.5:157.5) * pi/180;   % 8 个方向
    nPts  = 2*len + 1;
    intens = zeros(8, nPts);

    % 采样坐标
    t = linspace(-len, len, nPts);

    for k = 1:8
        th = theta(k);
        dx = cos(th);  dy = sin(th);
        x  = xc + t*dx;
        y  = yc + t*dy;
        intens(k,:) = interp2(1:W, 1:H, I, x, y, interpMethod);
    end

    % 可视化：原图叠加 8 条线
   % figHandle = figure; imshow(I, []); hold on;
   if iplot
        figHandle = figure('pos',[100,100,600,600]); 
        imagesc(I);axis equal; hold on;caxis([5e3,2e4]);
        colors = lines(8);
        for k = 1:8
            th = theta(k);
            dx = cos(th);  dy = sin(th);
            plot([xc-len*dx, xc+len*dx], [yc-len*dy, yc+len*dy], ...
                 'Color', colors(k,:), 'LineWidth', 1.5);
        end
        title('原图 + 8 条中心射线');
   end
  
end