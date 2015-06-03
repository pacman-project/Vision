% this is for 3D visualization of the objects represented in terms of
% the vocabulary

dataSetNumber = 3;
nClusters = 7;
layerID = 4;

% define folders configuration
commonRoot = 'D:/'; 
root = [commonRoot, 'LibHoP3D/'];
addPaths(root);

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

statistics1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
load(statistics1Layer);

displ3 = 6;
displ5 = 18;
displ7 = 52;


fileForVisualization{3} = [root, 'statistics/fileForVisualization_3_', dsN, '_', nCl, '.mat'];
fileForVisualization{4} = [root, 'statistics/fileForVisualization_4_', dsN, '_', nCl, '.mat'];
fileForVisualization{5} = [root, 'statistics/fileForVisualization_5_', dsN, '_', nCl, '.mat'];
fileForVisualization{6} = [root, 'statistics/fileForVisualization_6_', dsN, '_', nCl, '.mat'];
fileForVisualization{7} = [root, 'statistics/fileForVisualization_7_', dsN, '_', nCl, '.mat'];
fileForVisualization{8} = [root, 'statistics/fileForVisualization_8_', dsN, '_', nCl, '.mat'];

triple3OutDepth = [];
triple4OutDepth = [];
triple5OutDepth = [];
triple6OutDepth = [];
triple7OutDepth = [];
triple8OutDepth = [];

load(fileForVisualization{layerID});


% define a path to input files

[ depthPathDefault, depthPath ] = getPathToData(dataSetNumber, commonRoot, root);
fileListPrecomputed = false;
is_subset = false;
subset_len = 1;
subsetPercent = 0.01;

is_inhibition = true;

if is_inhibition
    % this are lists of files with detections:
    outRoot2 = [depthPath, '_layer2'];
    outRoot3 = [depthPath, '_layer3'];
    outRoot4 = [depthPath, '_layer4'];
    outRoot5 = [depthPath, '_layer5'];
    outRoot6 = [depthPath, '_layer6'];
    outRoot7 = [depthPath, '_layer7'];
    outRoot8 = [depthPath, '_layer8'];
    
end

depthPath = 'D:\Input Data\VladislavSTD\Vladislav_STD\depth_Inh';    
elPath = 'D:\Input Data\VladislavSTD\Vladislav_STD\depth_layer4_after_inh2';

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    [list_els, lenF] = extractFileList(fileListPrecomputed, elPath, '', is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subsetPercent);
    [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, subsetPercent);
end

[dxKernel, dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);

isY = true;
isX = true;
isTrim = false;


depthStep = thresh/4;
offsetD = 10;    
partColor = [0,0,1];
surfColor = [1,0,0];

isFlip = false;
if dataSetNumber == 2
    isFlip = true;
end

transparency = 0.3;

for i = 1:lenF
    I = imread(list_depth{i});
    I = I(:,:,1);
    marks = imread(list_els{i});
    
    [I, ~, ~, ~, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discRadius, isX, isY, false, ...
                            isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, [], [], [], []);
                        
    
%     imtool(I);
%     imtool(marks);
    
    maxRC = max([r,c]);
    [rows, cols] = find(marks>0);
    numEl = length(rows);
    depths = zeros(1,numEl);
    
    for j = 1:numEl
        % extract depths
        depths(j) = I(rows(j), cols(j));
    end
    
    
    minD = min(depths);
    maxD = max(depths);
    objectSizeDepth = ceil((maxD - minD)/depthStep);
    depthBottom = 3 * objectSizeDepth; % bottom surface for visualization
    depthRange =  7 * objectSizeDepth;
    
    f = figure;
    
    fieldSize = [maxRC, maxRC, depthRange];
    
    positions = [];
    elements = [];
    
    if layerID == 2  % simply make a list of elements
        elements = zeros(numEl,1);
        positions = zeros(numEl,3);
        for j = 1:numEl
            elements(j) =  marks(rows(j), cols(j));
            curDepth = floor((I(rows(j), cols(j)) - minD ) / depthStep) + depthBottom;
            positions(j,:) = [cols(j), rows(j), curDepth];
        end
    elseif layerID > 2
        for j = 1:numEl
            curDepth = floor((I(rows(j), cols(j)) - minD ) / depthStep) + depthBottom;
            elCenter = [cols(j), rows(j), curDepth];
            [poss, els] = partMeanReconstruction(layerID, marks(rows(j), cols(j)), elCenter, tripleOutDepth, displ3, displ5, displ7, nClusters);
            
            elements = [elements, els];
            positions = [positions; poss];
        end
    end

    [out] = surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, partColor, isFlip, transparency);
    
    isFIG = false;
    str_folder = 'D:\Visualized Reconstructions\png';
    A = exist(str_folder, 'dir');

    if A == 0
        mkdir(str_folder);
    end

    if ~ isFIG
        str = ['layer', num2str(layerID), '_', num2str(i), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');
    else
        str = ['layer', num2str(i), '.fig'];
        str1 = [str_folder, str];
        saveas(f, str1, 'fig');
    end
    

end



