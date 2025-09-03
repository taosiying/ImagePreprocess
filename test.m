% image = load("G:\Neutron Image Processing\imageReconEMML0608\recon0608\imagepreprocess\测试数据1&2&4\2-NISImageNew_p0_0.05_p2_0.25.mat").imageNIS;
% image = load("D:\works\imageprocess\NIS_data_20250530\p0_0.05_p2_0.0\Ideal\IdealImage_p0_0.05_p2_0.0.mat").imageNIS;
% image = load("D:\works\imageprocess\NIS_data_20250530\p0_0.05_p2_0.0\ColorNoise\ColorNoiseImage_5e2_p0_0.05_p2_0.0.mat").imageNIS;
% image = imread("F:\ImagePreProcess\datas\original_img\143.tif");
% image = double(image);
% r = 30;       % 滤波半径
% eps = 10; % 正则化参数
% denoised = guidedFilter(image, image, r, eps);

% denoised = savitzky_1d(image,51);
% denoised = waveletsdenoise(image);
% denoised = imguidedfilter(image);  
% figure();imagesc(image);title("original");axis equal tight;colorbar;
% figure();imagesc(denoisedeps);title("denoised");axis equal tight;colorbar;
% figure();imagesc(denoised - image);title("subtract");axis equal tight;colorbar;


% image1 = imread("D:\works\imageprocess\wzbGZDT\original_img\143.tif");
% [imageNISnew,G,indexIn]=imagepreprocess(image);
% imageNISnew1 = savitzky_1d(imageNISnew);
% [indexIn,imageNISnew]=imagepreprocess1(image);
% % imagesc(imageNISnew)
imageNIS = load("F:\ImagePreProcess\datas\TestExamples\TestExample1\NISImage_p0_0.05_p2_0.25_range0.1.mat").imageNIS;
rows = size(imageNIS, 1);
cols = size(imageNIS, 2);

noise_ratio = 0.16;
[u, v] = meshgrid(1:cols, 1:rows);
u = u - (cols + 1) / 2;
v = v - (rows + 1) / 2;
f = 1 ./ (1 + u.^2 + v.^2);
f = fftshift(f);  % 将滤波器的零频率移到中心
H = fft2(f); % 对滤波器进行傅里叶变换
noise = real(ifft2(H .* fft2(randn(rows, cols))));  % 生成色噪声
noise = noise * max(max(imageNIS))*noise_ratio; % 将噪声调整到合适的强度

imageNIS_noise = imageNIS + noise;
imageNISnew = waveletsdenoise(imageNIS_noise,5,7);
binary = seg_easily(imageNISnew);
xc,yc,a,b,theta = fit_ellipse_lsp(binary);
mask = get_ring(imageNISnew,0,601,601,a,b,theta);

% imageNISnew2 = savitzky_1d(imageNISnew,11);
% imageNISnew = waveletsdenoise(imageNISnew,11);
% imageNISnew = savitzky_1d(imageNISnew,31);
figure()
[M,N]=size(imageNISnew);
m = floor(M/2);
plot(imageNISnew(:,m));
hold on;
plot(imageNIS_noise(:,m));
hold on;
% plot(imageNISnew2(:,m),LineWidth=2);
% hold on;
% mask = zeros(M, N);
% mask(indexIn) = 2000;
plot(mask(:,m));
hold on;
% image = imresize(image,[M,N]);
% hold on;
% plot(imageNISnew1(:,m))
% hold on;
% imageideal = load("D:\works\imageprocess\NIS_data_20250530\p0_0.05_p2_0.0\Ideal\IdealImage_p0_0.05_p2_0.0.mat").imageNIS;
% plot(imageideal(:,m));
% hold on;
plot(imageNIS(:,m));
hold off;
