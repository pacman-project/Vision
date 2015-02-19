% this is to try arror_func


I = imread('concordaerial.png');
Igpu = gpuArray(I); 

Igray_gpu = arrayfun(@rgb2gray_custom,Igpu(:,:,1),Igpu (:,:,2),Igpu(:,:,3));

I_gpuresult = gather(Igray_gpu);

imtool(I_gpuresult);