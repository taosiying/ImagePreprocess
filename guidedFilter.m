function q = guidedFilter(I, p, r, eps)
% I: 引导图像（灰度，double，范围[0,1]）
% p: 输入图像（待滤波图像）
% r: 滤波半径
% eps: 正则化参数（防止除零）

% 计算局部均值
mean_I = imboxfilt(I, 2*r+1);
mean_p = imboxfilt(p, 2*r+1);

% 计算局部相关
corr_I = imboxfilt(I.*I, 2*r+1);
corr_Ip = imboxfilt(I.*p, 2*r+1);

% 计算方差和协方差
var_I = corr_I - mean_I .* mean_I;
cov_Ip = corr_Ip - mean_I .* mean_p;

% 计算线性系数 a 和 b
a = cov_Ip ./ (var_I + eps);
b = mean_p - a .* mean_I;

% 计算均值 a 和 b
mean_a = imboxfilt(a, 2*r+1);
mean_b = imboxfilt(b, 2*r+1);

% 输出滤波结果
q = mean_a .* I + mean_b;
end
