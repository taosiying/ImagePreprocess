imgRaw = imread("D:\works\nispreprocess\original_img\145.tif");
background = imread("D:\works\nispreprocess\background.tif");
smooth_window = 51;sample_ratio = 0.01;dr = 200;save_name = '145';
imageNIS = nispreprocess(imgRaw,background,smooth_window,sample_ratio,dr,save_name);
RadialProfile_cmp(imgRaw,imageNIS,2000,2000,0,0,"处理前后对比",false);