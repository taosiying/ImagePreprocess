%%
% 减去本底主函数:首先分别计算四角暗区和半影区，再对背景区拟合一次曲面1，减去曲面1后对中心平顶区拟合曲面2，减去曲面2后得到最终结果
% I:输入图像
% sample_ratio:拟合数据对输入图像的采样比例，默认为0.1,实际图像可调整为0.01
% dr：划分半影区向内外扩充的半径大小，默认为100，实际图像需要调整为200，防止半影区的趋势影响拟合曲面趋势
% 输出I_subed2:去本底后曲面，surface_sum:拟合后曲面
% 使用示例：imgBkg = sub_surface(imgSmooth,0.01,15,200);
%%

function [I_subed2,surface_sum] = sub_surface(I,sample_ratio,dr)
if nargin < 2
    sample_ratio = 0.1;
end
if nargin < 3
     dr = 100;
end
I1 = savitzky_1d(I,101);
I1 = savitzky_1d(I1,101);
binary = seg_easily(I1);
[xc,yc,a,b,theta] = fit_ellipse_lsp( binary );
[~,outside_mask,inside_mask] = get_ring(I,dr, xc, yc, a, b, theta);
mask1 = zeros(size(I1));
image_mask = (I - min(I(:)))/(max(I(:))-min(I(:)));
mask1(image_mask<0.2) = 255;
mask1 = uint8(mask1);
se = strel('disk', 31);      
mask1 = imdilate(mask1, se);
mask1 = uint8(imbinarize(mask1))*255;
mask1 = bitcmp(mask1);
mask =bitand(mask1,bitcmp(outside_mask));
surface = fit_gauss(I1,mask,5,sample_ratio,100,2,false);
I_subed = I - surface;
I_subed = I_subed- mean(I_subed(mask>0));
surface2 = fit_gauss(I_subed,inside_mask,5,sample_ratio,100,2,true);

I_subed2 = I_subed - surface2 + min(surface2);
I_subed2 = I_subed2 + mean(I1(inside_mask>0)) - mean(I_subed2(inside_mask>0));
I_subed2 = I_subed2 - mean(I_subed2(mask>0));
I_subed2(mask1==0) = 0;
surface_sum = I - I_subed2;
figure();
subplot(2,3,1);
imagesc(I1);axis image; colorbar;title("平滑后曲面");
subplot(2,3,2);
imagesc(surface);axis image; colorbar;title("背景区拟合曲面");
subplot(2,3,3);
imagesc(surface2);axis image; colorbar;title("平顶区拟合曲面");
subplot(2,3,4);
imagesc(I);axis image; colorbar;title("原图");
subplot(2,3,5);
imagesc(I_subed);axis image; colorbar;title("减去背景区拟合曲面后");
subplot(2,3,6);
imagesc(I_subed2);axis image; colorbar;title("减去平顶区拟合曲面后");

end

