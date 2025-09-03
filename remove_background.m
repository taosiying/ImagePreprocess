%%
% 平场校正函数，输入图像img和校正图像background，
% 首先对background进行归一化，对background小于0.3的值置为1避免极值影响，然后将img与background相除，再除以系数保持均值不变，得到校正结果
% 使用示例：
% background = imread("background.tif");
% imgFlat = remove_background(imgRaw,background);

%%
function img = remove_background(img,background)
    % Remove background using a reference image
    min_img = min(img(:));
    img = img - min_img;
    mean_inside_original = mean(img(:));
    disp(mean_inside_original);
    [height, width] = size(img); 
    background = imresize(background, [height, width]);
    background = double(background);
    background = (background - min(background(:))) / (max(background(:)) - min(background(:)));
    background(background<0.3) = 1;%避免极值影响
    img = double(img) ./ background;
    mean_inside_after = mean(img(:));
    img = img ./ mean_inside_after .* mean_inside_original; %保持均值不变
    img = img + min_img;
end