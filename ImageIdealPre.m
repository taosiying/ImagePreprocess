%% 主脚本：中子成像处理流程
%% 理想图像处理脚本：包含加本底，加噪声，去噪，去本底
clc; close all; clear;

%% 1. 读图
Dir= 'E:\ImagePreProcess\datas\';
[filename, filepath] = uigetfile(fullfile(Dir,'*.mat'),'选择.mat文件'); % 交互选 .mat、.tif
data = load(fullfile(filepath, filename));
imgRaw   = double(load(fullfile(filepath, filename)).imageNIS);
imgRaw = imgRaw - min(imgRaw(:));
[path, name, ext] = fileparts(fullfile(filepath, filename));
dirname=[pwd,'\',name];
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
A_background = max(imgRaw(:)) * 0.3;
A_gaussian = A_background;
num_gaussian = 3;
img_size = size(imgRaw,1);
background = generate_background(img_size,A_background,A_gaussian,num_gaussian);
img_background = background + imgRaw;
noise_ratio = 0.05;
[rows,cols] = size(img_background);
[u, v] = meshgrid(1:cols, 1:rows);
u = u - (cols + 1) / 2;
v = v - (rows + 1) / 2;
f = 1 ./ (1 + u.^2 + v.^2);
f = fftshift(f);  % 将滤波器的零频率移到中心
H = fft2(f); % 对滤波器进行傅里叶变换
noise = real(ifft2(H .* fft2(randn(rows, cols))));  % 生成色噪声
noise = noise * max(max(img_background))*noise_ratio; % 将噪声调整到合适的强度
img_noise = img_background + noise;
figure; imagesc(imgRaw); axis image;%caxis([5e3,2e4]);
figure; imagesc(img_background); 
figure; imagesc(img_noise);
xc = 601;yc = 601;
%% 2. 小波去噪
% [thr, sorh, keepapp] = ddencmp('den', 'wv', imgRaw);
% imgDen = wdencmp('gbl', imgRaw, 'sym4', 4, thr, sorh, keepapp);
imgDen = waveletsdenoise(img_noise);
    %检测
    string='小波去噪';
    RadialProfile_cmp(img_noise,imgDen,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);

%% 4. 匀滑
% imgFlat(~isfinite(imgFlat)) = 0;
imgSmooth = savitzky_1d(imgDen, 13);
    %检测
    string='图像匀滑';
    RadialProfile_cmp(imgDen,imgSmooth,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
%% 5. 去本底（圆外置 0 → 圆内置 0 → 曲面拟合 → 减本底）
    [imgBkg,surface_pre] = sub_surface(imgSmooth,0.1,100);
    string='去本底';
    RadialProfile_cmp(imgSmooth,imgBkg,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
    string='和理想图像对比';
    RadialProfile_cmp(imgRaw,imgBkg,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
%% 6. 二次手工圆掩膜（外围置 0）
% disp('边缘背景归零：交互圆外置 0');
% [imgCrop, ~, ~] = circleMaskManual(imgBkg,[5e3,2e4]);  % 交互裁剪
% imgCrop(~isfinite(imgCrop)) = 0;
%         %检测
%          string='二次掩膜';
%         RadialProfile_cmp(imgBkg,imgCrop,xc,yc,0,0,string);
%         saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
%         drawnow;pause(3);
imgCrop = imgBkg;
%% 7. 坐标单位重构（mm）
px = 36.9e-3;                          % 像素尺寸 mm
xPix = (1:size(imgCrop,2)) - (size(imgCrop,2)+1)/2;
yPix = (1:size(imgCrop,1)) - (size(imgCrop,1)+1)/2;
xMM  = xPix * px;
yMM  = yPix * px;

imageNIS  = imgCrop;
imagePlane = [min(xMM), max(xMM), min(yMM), max(yMM), mean(diff(xMM))];

%% 8. 保存结果
% save([dirname, '\',name,'.mat'], 'imageNIS', 'imagePlane', 'xMM', 'yMM');
% save([dirname, '\',name,'_proData.mat'], 'imgRaw', 'imgDen','imgSmooth','imgBkg','imageNIS');
Image_out.imageNIS = imageNIS;
Image_out.xList = data.xList;
Image_out.yList = data.yList;
save([dirname, '\',name,'.mat'],"Image_out");
%% 9. 可视化

figure('Position', [100 100 1500 900]);
subplot(2,3,1);imagesc(imgRaw); axis equal;title('原始数据：imgRaw');
subplot(2,3,2);imagesc(img_background);axis equal;title('加入背景：imgbackground');
subplot(2,3,3);imagesc(img_noise);axis equal;title('加入噪声：imgnoise');
subplot(2,3,4);imagesc(imgDen); axis equal; axis equal;title('去噪：imgDen');
subplot(2,3,5);imagesc(imgSmooth); axis equal; axis equal;title('平滑：imgSmooth');
subplot(2,3,6);imagesc(imgBkg); axis equal; axis equal;title('去背景：imgBkg');
diff_image =abs(imgRaw - imgBkg);
string='数据处理流程';
saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
figure();imagesc(diff_image);colorbar;title("与理想图像差值");



figure; imagesc(xMM, yMM, imageNIS); axis equal tight;
xlabel('X (mm)'); ylabel('Y (mm)');
title('Measured Image (mm scale)');
string='最终结果：imageNIS';
saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);

% SS_residual = sum(sum( (diff_image).^2 ));
% SS_total = sum(sum( (img_noise- mean(img_noise)).^2 ));
% 
% % 5. 计算决定性系数 R²
% R2 = 1 - (SS_residual / SS_total);
% disp(R2)
