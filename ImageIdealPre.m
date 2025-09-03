%% ���ű������ӳ���������
%% ����ͼ����ű��������ӱ��ף���������ȥ�룬ȥ����
clc; close all; clear;

%% 1. ��ͼ
Dir= 'E:\ImagePreProcess\datas\';
[filename, filepath] = uigetfile(fullfile(Dir,'*.mat'),'ѡ��.mat�ļ�'); % ����ѡ .mat��.tif
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
f = fftshift(f);  % ���˲�������Ƶ���Ƶ�����
H = fft2(f); % ���˲������и���Ҷ�任
noise = real(ifft2(H .* fft2(randn(rows, cols))));  % ����ɫ����
noise = noise * max(max(img_background))*noise_ratio; % ���������������ʵ�ǿ��
img_noise = img_background + noise;
figure; imagesc(imgRaw); axis image;%caxis([5e3,2e4]);
figure; imagesc(img_background); 
figure; imagesc(img_noise);
xc = 601;yc = 601;
%% 2. С��ȥ��
% [thr, sorh, keepapp] = ddencmp('den', 'wv', imgRaw);
% imgDen = wdencmp('gbl', imgRaw, 'sym4', 4, thr, sorh, keepapp);
imgDen = waveletsdenoise(img_noise);
    %���
    string='С��ȥ��';
    RadialProfile_cmp(img_noise,imgDen,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);

%% 4. �Ȼ�
% imgFlat(~isfinite(imgFlat)) = 0;
imgSmooth = savitzky_1d(imgDen, 13);
    %���
    string='ͼ���Ȼ�';
    RadialProfile_cmp(imgDen,imgSmooth,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
%% 5. ȥ���ף�Բ���� 0 �� Բ���� 0 �� ������� �� �����ף�
    [imgBkg,surface_pre] = sub_surface(imgSmooth,0.1,100);
    string='ȥ����';
    RadialProfile_cmp(imgSmooth,imgBkg,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
    string='������ͼ��Ա�';
    RadialProfile_cmp(imgRaw,imgBkg,xc,yc,0,0,string);
    saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
    drawnow;pause(3);
%% 6. �����ֹ�Բ��Ĥ����Χ�� 0��
% disp('��Ե�������㣺����Բ���� 0');
% [imgCrop, ~, ~] = circleMaskManual(imgBkg,[5e3,2e4]);  % �����ü�
% imgCrop(~isfinite(imgCrop)) = 0;
%         %���
%          string='������Ĥ';
%         RadialProfile_cmp(imgBkg,imgCrop,xc,yc,0,0,string);
%         saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);
%         drawnow;pause(3);
imgCrop = imgBkg;
%% 7. ���굥λ�ع���mm��
px = 36.9e-3;                          % ���سߴ� mm
xPix = (1:size(imgCrop,2)) - (size(imgCrop,2)+1)/2;
yPix = (1:size(imgCrop,1)) - (size(imgCrop,1)+1)/2;
xMM  = xPix * px;
yMM  = yPix * px;

imageNIS  = imgCrop;
imagePlane = [min(xMM), max(xMM), min(yMM), max(yMM), mean(diff(xMM))];

%% 8. ������
% save([dirname, '\',name,'.mat'], 'imageNIS', 'imagePlane', 'xMM', 'yMM');
% save([dirname, '\',name,'_proData.mat'], 'imgRaw', 'imgDen','imgSmooth','imgBkg','imageNIS');
Image_out.imageNIS = imageNIS;
Image_out.xList = data.xList;
Image_out.yList = data.yList;
save([dirname, '\',name,'.mat'],"Image_out");
%% 9. ���ӻ�

figure('Position', [100 100 1500 900]);
subplot(2,3,1);imagesc(imgRaw); axis equal;title('ԭʼ���ݣ�imgRaw');
subplot(2,3,2);imagesc(img_background);axis equal;title('���뱳����imgbackground');
subplot(2,3,3);imagesc(img_noise);axis equal;title('����������imgnoise');
subplot(2,3,4);imagesc(imgDen); axis equal; axis equal;title('ȥ�룺imgDen');
subplot(2,3,5);imagesc(imgSmooth); axis equal; axis equal;title('ƽ����imgSmooth');
subplot(2,3,6);imagesc(imgBkg); axis equal; axis equal;title('ȥ������imgBkg');
diff_image =abs(imgRaw - imgBkg);
string='���ݴ�������';
saveas(gcf, [dirname,'\',string,'.png']);%saveas(gcf, [dirname,'\',string,'.fig']);
figure();imagesc(diff_image);colorbar;title("������ͼ���ֵ");



figure; imagesc(xMM, yMM, imageNIS); axis equal tight;
xlabel('X (mm)'); ylabel('Y (mm)');
title('Measured Image (mm scale)');
string='���ս����imageNIS';
saveas(gcf, [dirname,'\',string,'.png']);saveas(gcf, [dirname,'\',string,'.fig']);

% SS_residual = sum(sum( (diff_image).^2 ));
% SS_total = sum(sum( (img_noise- mean(img_noise)).^2 ));
% 
% % 5. ���������ϵ�� R�0�5
% R2 = 1 - (SS_residual / SS_total);
% disp(R2)
