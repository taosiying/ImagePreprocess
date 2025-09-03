%% 主脚本：中子成像处理流程
%% 实际图像处理脚本，包含平场校正、去噪、平滑、去本底
clc; close all; clear;

%% 1. 读图
Dir= 'D:\works\nispreprocess\original_img';
[filename, filepath] = uigetfile(fullfile(Dir,'*.tif'),'选择.tif文件'); % 交互选 .mat、.tif

imgRaw   = double(imread(fullfile(filepath, filename)));
[path, name, ext] = fileparts(fullfile(filepath, filename));
dirname=[pwd,'\',name];
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
% binary = seg_easily(imgRaw);
% [xc,yc,a,b,theta] = fit_ellipse_lsp(binary);
figure; imagesc(imgRaw); axis image;%caxis([5e3,2e4]);
[xc, yc] = ginput(1);   % 点一次
%% 2. 平场修正
background = imread("background.tif");
imgFlat = remove_background(imgRaw,background);
% imgFlat = imgMed ./ imFlat;
    %检测
        string='平场修正';
        RadialProfile_cmp(imgRaw,imgFlat,xc,yc,0,0,string);
        saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
        drawnow;pause(3);
%% 3. 去噪
imgDen = waveletsdenoise(imgFlat);
    %检测
    string='去噪';
    RadialProfile_cmp(imgFlat,imgDen,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);

%% 4. 匀滑
imgDen(~isfinite(imgFlat)) = 0;
% imgSmooth = Fn_smooth2D(imgFlat, 5);
imgSmooth = savitzky_1d(imgDen,51);
    %检测
    string='图像匀滑';
    RadialProfile_cmp(imgDen,imgSmooth,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
%% 5. 去本底（圆外置 0 → 圆内置 0 → 曲面拟合 → 减本底）
imgBkg = sub_surface(imgSmooth,0.01,200);
    %检测
        string='去本底';
        RadialProfile_cmp(imgSmooth,imgBkg,xc,yc,0,0,string);
        saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
        drawnow;pause(3);
%% 6. 二次手工圆掩膜（外围置 0）
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
save([dirname, '\',name,'.mat'], 'imageNIS', 'imagePlane', 'xMM', 'yMM');
save([dirname, '\',name,'_proData.mat'], 'imgRaw', 'imgDen',...
    'imgFlat','imgSmooth','imgBkg','imgCrop','imageNIS');
%% 9. 可视化

% 绘制 8 幅子图
figure('Position', [100 100 1500 900]);
subplot(2,3,1);imagesc(imgRaw);caxis([5e3 2e4]); axis equal;title('原始数据：imgRaw');colorbar;
subplot(2,3,2);imagesc(imgFlat);caxis([5e3 2e4]); axis equal; axis equal;title('平场修正：imgFlat');colorbar;
subplot(2,3,3);imagesc(imgDen); axis equal;title('去噪：imgDen');colorbar;
subplot(2,3,4);imagesc(imgSmooth); axis equal; axis equal;title('高斯平滑：imgSmooth');colorbar;
subplot(2,3,5);imagesc(imgBkg); axis equal;title('去背景：imgBkg');colorbar;
subplot(2,3,6);imagesc(imageNIS);  axis equal;title('最终结果：imageNIS');colorbar;
string='数据处理流程';
saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
% 
% figure; imagesc(xMM, yMM, imageNIS); axis equal tight;
% xlabel('X (mm)'); ylabel('Y (mm)');
% title('Measured Image (mm scale)');
% caxis([0 8000]);
% string='最终结果：imageNIS';
% saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
% 
% figure; imagesc(imgRaw);caxis([5e3 2e4]); axis equal;
% string='原始数据：imgRaw';title(string);
% saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
% 
% figure; imagesc(imgDen);caxis([5e3 2e4]); axis equal;
% string='小波去噪：imgDen';title(string);
% saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
% 
% 
% figure; imagesc(imgFlat);caxis([5e3 2e4]); axis equal;
% string='平场修正：imgFlat';title(string);
% saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
% 
% figure; imagesc(imgSmooth);caxis([5e3 2e4]); axis equal;
% string='平滑：imgSmooth';title(string);
% saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
% 
% figure; imagesc(imgBkg);caxis([5e3 2e4]); axis equal;
% string='去背景：imgBkg';title(string);
% saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
% 
% figure; imagesc(imgCrop);caxis([5e3 2e4]); axis equal;
% string='边缘背景归零：imgCrop';title(string);
% saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
