% this is for 3D visualization of the objects represented in terms of
% vocabulary

% define folders configuration
root = 'D:\3D\LibHoP3D\';
addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'Optimization']);
addpath([root,'Learning']);
addpath([root,'visualiaztion']);

dataRoot = 'D:\3D\Elements2Layer\';
dataEl = 'Elements1';
dataTI = 'TrimedImages1';

lineDefault = [];

fileListPrecomputed = false; 
is_subset = false; % whether we shall use all files for learning
subset_len = 1000; % how much shall we use
nClusters = 9;

[list_marks, lenF] = extractFileList(fileListPrecomputed, [dataRoot, dataEl], lineDefault, is_subset, subset_len);
[list_TI, lenF] = extractFileList(fileListPrecomputed, [dataRoot, dataTI], lineDefault, is_subset, subset_len);

struct = load([root, 'settings/firstLayer']); % read the first layer coverage
cluster1Centres = struct.cluster1Centres;
cluster1Lengths = struct.cluster1Lengths;
thresh = struct.thresh;

    

for i = 1:1 % lenF
    I = imread(list_TI{i});
    marks = imread(list_marks{i});
    
    imtool(I);
    imtool(marks);
    
    % define the depth range as 10 percent quantile
    a = I(:);
    minDepth = quantile(a, 0.1);
    maxDepth = max(a);
    [r,c] = size(I);
    depthStep = thresh/3;
    
    %fieldSize = [c,r, round((maxdepth - minDepth)/depthStep)];
    
    maxRC = max([r,c]);
    
    fieldSize = [maxRC, maxRC, round((maxDepth - minDepth)/depthStep)];

    [rows, cols] = find(marks>0);
    numEl = length(rows);
    elements = zeros(1,numEl);
    positions = zeros(numEl, 3);
    
    for j = 1:numEl
        elements(j) = marks(rows(j), cols(j));
        curDepth = floor(I(rows(j), cols(j)) / depthStep);
        if curDepth < minDepth 
            curDepth = minDepth;
        end
        if curDepth > maxDepth 
            curDepth = maxDepth;
        end
        positions(j,:) = [cols(j),rows(j),curDepth];
    end
    
    [out] = surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
end



