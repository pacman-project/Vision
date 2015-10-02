clc;
close all;
clear all;

theta = [0:15:360];
inImg = imread('example.png');
dim = ndims(inImg);
if(dim == 3)
    %Input is a color image
    inImg = rgb2gray(inImg);
end
I = inImg;
tic
for i = [1:length(theta)]
   J1(:,:,i) = steerGaussFilterOrder1(I,theta(i),4, true);
end
toc

tic
for i = [1:length(theta)]
   J2(:,:,i) = steerGaussFilterOrder2(I,theta(i),3,true);
end
toc