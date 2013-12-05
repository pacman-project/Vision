% this is a script for work with Mete and Mirela
% work with

dataRoot = 'D:/3D/Input Data/washington_dataset/';
dataFolder = 'coffee_mug_';

root = 'D:/3D/LibHoP3D/';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'Optimization']);
addpath([root,'Learning']);
addpath([root,'washingtonSettings']);

list_depth = {};

% creating file list
for i = 1:8
    curPath = [dataRoot, dataFolder, num2str(i)];
    cur_list_depth = load_filelist(curPath);
    list_depth = [list_depth, cur_list_depth];
    i
end

ss = load('settings/sigma.mat');
sigma = ss.sigma; % this is a sigma for gaussian derivative 
sigmaKernelSize = ss.sigmaKernelSize; % size of the kernel for gaussian
dxKernel = load('dxKernel.mat');
dxKernel = dxKernel.dxKernel;
nClusters = 9;

[cluster1Centres, cluster1Lengths, thresh] = learnFirstLayer(list_depth, sigma, sigmaKernelSize, dxKernel, nClusters);



