function [imgMask, xC, yC, r] = circleMaskManual(img, Icolor,invert)
% invert=false: Բ����1��invert=true: Բ����1
if nargin < 3, invert = false; end
figure; imagesc(img);caxis(Icolor); 
axis equal;hold on;
title('1. ���Բ��'); [xC, yC] = ginput(1); plot(xC, yC, 'ro');
title('2. ���Բ��'); [xR, yR] = ginput(1); r = norm([xR-xC, yR-yC]);
[X, Y] = meshgrid(1:size(img,2), 1:size(img,1));
mask = ((X - xC).^2 + (Y - yC).^2) <= r.^2;
if invert, mask = ~mask; end
imgMask = img .* mask;

imagesc(imgMask);caxis(Icolor); 
axis equal;
%close(gcf);
end