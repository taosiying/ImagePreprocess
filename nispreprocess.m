% %% 主脚本：中子成像处理函数
% imgRaw：输入图像矩阵
% backround：平场校正矩阵，当不输入时，不进行平场校正处理
% smooth_window: 图像平滑窗口大小，小于3时不进行图像平滑
% sample_ratio:采样频率，默认0.01
% dr:半影区向内外扩张尺度，默认为100，实际图像建议设置为200
% savename:保存文件名，将保存为再当前目录下的savename.mat

function [imageNIS] = nispreprocess(imgRaw,background,smooth_window,sample_ratio,dr,save_name)
imgRaw = double(imgRaw);
if nargin < 6
    save_name = 'ImageNIS';
end
if nargin < 5
    dr = 100;
end
if nargin < 4
    sample_ratio = 0.01;
end
if nargin < 3
    smooth_window = 51;
end
if nargin < 2
    imgFlat = imgRaw;
else
    %% 平场校正
    disp("正在平场校正");
    imgFlat = remove_background(imgRaw,background);
end
    %% 3. 去噪
    disp("正在去噪");
    figure();imagesc(imgFlat);title("平场校正后图像");
    imgDen = waveletsdenoise(imgFlat);
    figure();imagesc(imgDen);title("去噪后图像");
    %% 4. 图像匀滑
    disp("正在平滑图像");
    if smooth_window>1
        imgSmooth = savitzky_1d(imgDen,smooth_window);
    else
        imgSmooth = imgDen;
    end
    %% 5. 去本底
    disp("正在去除本底");
    imgBkg = sub_surface(imgSmooth,sample_ratio,dr);
    imageNIS  = imgBkg;
    %% 保存
    dirname = pwd;
    save([dirname, '\',save_name,'.mat'], 'imageNIS');
    %% 9. 可视化
end