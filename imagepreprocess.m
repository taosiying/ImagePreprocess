%% 图像预处理
function [imageNISnew,G,indexIn]=imagepreprocess(imageNIS)
%%坐标变换参数
[X, Y] = meshgrid(1:size(imageNIS,1), 1:size(imageNIS,2));
xListTmp = X(1,:) - X(end)/2;
xListTmp = xListTmp*9*150/40/1e3;
yListTmp = xListTmp;
xList = fix(xListTmp(1)):0.04:fix(xListTmp(end));
yList = xList;
ispacevar = 0;
idetector = 1;
idiff = 0;
offShift = 0;
tiltAng = 0;

Lmfp=28.76;%mm

bin = 1;
imageP = [-60 60 -60 60 0.1*bin];
sourceP = [-0.06 0.06 -0.06 0.06 0.004]*bin; %mm
%%
center0X=0;
center0Y=0;

drawnow;

uDataOri=[xListTmp(1) xListTmp(end)];
vDataOri=[yListTmp(1) yListTmp(end)];
xDataOri=uDataOri;
yDataOri=vDataOri;
xyScaleOri=(xDataOri(2)-xDataOri(1))/(size(imageNIS, 1)-1);

uData=[xList(1) xList(end)];
vData=[yList(1) yList(end)];
xData=imageP(1:2);
yData=imageP(3:4);
xyScale=imageP(5);
rotMatZ = rotz(0);
xListNew = xData(1):xyScale:xData(2);
yListNew = yData(1):xyScale:yData(2);
[xListNewMesh, yListNewMesh] = meshgrid(xListNew, yListNew);
rListNew = sqrt(xListNewMesh.^2 + yListNewMesh.^2);
%%
figure;imagesc(xListTmp, yListTmp, imageNIS);axis equal tight;title('实验原始编码图像');colorbar;
%% 平场校正
corrected_image = remove_background(imageNIS);
denoised_image = waveletsdenoise(corrected_image);
figure();imagesc(xListTmp, yListTmp,denoised_image);title("denoised");axis equal tight;colorbar;

figure();imagesc(xListTmp, yListTmp,corrected_image);title("corrected_image");axis equal tight;colorbar;
%% 二值化
binary = seg_easily(denoised_image);
% figure();imshow(binary);title(binary);axis equal tight;
%% 拟合椭圆
[xc,yc,a,b,theta] = fit_ellipse_lsp( binary );
center0X = (xc - size(imageNIS, 2)/2) /2048*mean(abs(uDataOri));
center0Y = (yc - size(imageNIS, 1)/2) /2048*mean(abs(uDataOri));
% fprintf("圆心：(%f,%f),半径：(%f,%f)",xc,yc,a,b)
%% 图像平滑
smoothed_image =  savitzky_1d(denoised_image, 51);
figure();plot(denoised_image(:,size(imageNIS, 1)/2));
hold on;plot(smoothed_image(:,size(imageNIS, 1)/2));
% figure();imagesc(xListTmp, yListTmp,smoothed_image);title("smoothed_image");axis equal tight;colorbar;
%% 减去本底
% suf_background_image = sub_background(smoothed_image, xc, yc, a, b, theta);
suf_background_image = sub_surface(smoothed_image);
dr = 50;
 T=maketform('affine',rotMatZ);
 [imageNIS,~,~]=imtransform(suf_background_image,T,'bicubic','UData',uDataOri - center0X,'VData',...
    vDataOri - center0Y,'XData',xDataOri,'YData',yDataOri,'XYScale',xyScaleOri,'FillValues',0);
figure,imagesc(xListTmp, yListTmp, imageNIS);axis equal tight;title('suf background centered_mage');

[XmeshTmp, YmeshTmp] = meshgrid(xListTmp, yListTmp);
[Xmesh, Ymesh] = meshgrid(xList, yList);
imageNIS2 = interp2(XmeshTmp, YmeshTmp, imageNIS, Xmesh, Ymesh);
 [imageNISnew,~,~]=imtransform(imageNIS2,T,'bicubic','UData',uData,'VData',...
    uData,'XData',xData,'YData',yData,'XYScale',xyScale,'FillValues',0);
binary = seg_easily(imageNISnew);
% figure();imshow(binary);
[xc,yc,a,b,theta] = fit_ellipse_lsp(binary);
dr = 10;
%% 计算半影区
ring_mask = get_ring(imageNISnew,dr,xc,yc,a,b,theta);
figure();imshow(ring_mask);

figure,imagesc(xData, yData, imageNISnew);axis equal tight;
hold on;

%% 绘制 ring_mask 的边缘
boundaries = bwboundaries(ring_mask);
for k = 1:length(boundaries)
    boundary = boundaries{k};
    boundary_x = (boundary(:, 2) - size(imageNISnew, 2)/2) * xyScale;
    boundary_y = (boundary(:, 1) - size(imageNISnew, 1)/2) * xyScale;
    plot(boundary_x, boundary_y, 'y-', 'LineWidth', 1.);
end
title('划分ROI区间后的实验编码图像');
drawnow;
%% G计算
r0NIS=(a+b)/2;
G=(r0NIS-17.83)/6.761;
fprintf("G:%f",G)
indexIn = find(ring_mask>0);
end
