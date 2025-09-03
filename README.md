# ImagePreprocess
中子图像预处理算法

## 主函数：`nispreprocess.m`

### 输入
- `imgRaw`：输入图像矩阵  
- `background`：平场校正矩阵（可选，不输入时不进行平场校正处理）  
- `smooth_window`：图像平滑窗口大小，小于 3 时不进行平滑  
- `sample_ratio`：采样频率，默认 0.01  
- `dr`：半影区向内外扩张尺度，默认 100，实际图像建议设置为 200  
- `savename`：保存文件名，将保存为当前目录下的 `savename.mat`

### 输出
- `imageNIS`：处理后的图像

### 使用示例
```matlab
imgRaw = imread("D:\works\nispreprocess\original_img\145.tif");
background = imread("D:\works\nispreprocess\background.tif");
smooth_window = 51;
sample_ratio = 0.01;
dr = 200;
save_name = '145';
imageNIS = nispreprocess(imgRaw, background, smooth_window, sample_ratio, dr, save_name);
```
## 对实际图像处理主脚本：`ImageIdealPre.m`
-实际图像处理脚本，包含平场校正、去噪、平滑、去本底
-对每一步处理后输出并通过灰度曲线验证

## 对理想图像验证主脚本：ImageIdealPre.m
包含加本底，加噪声，去噪，去本底, 对每一步处理后输出并通过灰度曲线验证
图像平场校正函数：remove_background.m
平场校正函数，首先对background进行归一化，对background小于0.3的值置为1避免极值影响，然后将img与background相除，再除以系数保持均值不变，得到校正结果
### 输入：
图像img和校正图像background
### 输出：
校正后图像
### 使用示例：
```matlab
background = imread("background.tif");
imgFlat = remove_background(imgRaw,background);
图像去噪函数：waveletsdenoise.m
```
## 基于小波变换的去噪函数：waveletsdenoise.m
将图像分解为四个频带，分别为低频AA自带，HL子带，LH自带，HH子带，在AA子带保持不变，在其它子带分别用sobel算子计算边缘，在边缘区保持不变，其它区域通过软阈值变换，去除噪声。
### 输入：
noisy_img:带噪声图像，window_size1:对低频自带平滑窗口(已禁用，保持低频自带不变)，window_size2:对整体图像平滑窗口大小
### 输出：
校正后图像
### 使用示例：
```matlab
imgDen = waveletsdenoise(imgFlat，5，7);
```
##图像平滑函数：savitzky_1d.m
### 输入： 
img：输入图像，window_size:平滑窗口大小
### 使用示例：
imgSmooth = savitzky_1d(imgDen,51);（实际图像）
imgSmooth = savitzky_1d(imgDen,31);（仿真图像）
## 曲面拟合和去背景函数：sub_surface.m
减去本底主函数:首先分别计算四角暗区和半影区，再对背景区拟合一次曲面1，减去曲面1后对中心平顶区拟合曲面2，减去曲面2后得到最终结果
### 输入： 
I: 输入图像
sample_ratio: 拟合数据对输入图像的采样比例，默认为0.1,实际图像可调整为0.01
dr：划分半影区向内外扩充的半径大小，默认为100，实际图像需要调整为200，防止半影区的趋势影响拟合曲面趋势
### 输出：
I_subed2:去本底后曲面，surface_sum:拟合后曲面
### 使用示例：
imgBkg = sub_surface(imgSmooth,0.01,15,200);
## 拟合曲面函数：fit_gauss.m
对图像image的mask区域的点按照sample_ratio比例采样，并通过最小二乘法，拟合为num_gaussians个高斯函数组成的曲面
### 输入：
-`image`：输入图像数据
-`mask`：被拟合区域二值图像
-num_gaussians:高斯函数个数
-sample_ratio:采样频率
-max_iterations:最大迭代次数
-max_sigma：最大方差比例，限制方差为max_sigma*size(image,1)
-constrain_flat:是否限制mask以外区域尽量平坦
### 使用示例：
```matlab
surface = fit_gauss(image,mask,5,sample_ratio,100,2,false);
```
## 生成随机高斯曲面本经的函数：generate_background.m
### 输入:
img_size:图像尺寸，A_background：背景圆幅度，A_gaussian:高斯曲面幅度，num_gaussian:高斯波包个数
### 输出：
与img_size相同的图像矩阵
### 使用示例：
```matlab
A_background = max(imgRaw(:)) * 0.3;
A_gaussian = A_background;
num_gaussian = 3;
img_size = size(imgRaw,1);
background = generate_background(img_size,A_background,A_gaussian,num_gaussian);
```

