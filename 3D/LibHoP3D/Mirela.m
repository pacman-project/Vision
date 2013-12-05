% this is a script for work with Mete and Mirela
% work with

dataRoot = 'F:/washington_dataset/';
dataFolder = 'coffee_mug_';

root = 'D:/3D/LibHoP3D/';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'Optimization']);
addpath([root,'Learning']);
addpath([root,'washingtonSettings']);

list_depth = {};

% creating file list
for i = 1:7
    curPath = [dataRoot, dataFolder, num2str(i)];
    cur_list_depth = load_filelist(curPath);
    list_depth = [list_depth, cur_list_depth];
    i
end

outFileList = {};
lenF = length(list_depth);
% inds = 1:3:lenF;
% list_depth{inds} = [];
% lenF = length(list_depth);

% apply mask to the depth crop images, remove crops and masks from the filelist
ss = num2str(lenF);

for i = 1:3:lenF
    % take all three corresponding files
    I = imread(list_depth{i+1});
    mask = imread(list_depth{i+2});
    mask = uint16(mask);
    I = I.*mask;
    str = list_depth{i+2};
    str = ['D:/3D/Input Data/', str(4:end)];
    imwrite(I, str, 'png');
    
    if mod(i, 10) == 0
        strOut = [num2str(i),'/',ss];
        disp(strOut);
    end
end

a = 2;
% [cluster1Centres, cluster1Lengths, thresh] = learnFirstLayer(list_depth, sigma, sigmaKernelSize, dxKernel, nClusters)
% [allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel)


