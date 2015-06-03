% this code is to create a MULTISCALE histogram of compositional parts

function [] = histogramOfParts_best_Washington(nNClusters, numChannels)

% pyramid_levels = [1.0, 0,66];
% pyramid_num_levels = length(pyramid_levels);

nbinsS1 = [1,1];
nbinsS2 = [2,2];
nbinsS3 = [3,3];

dataSetNumber = 2;
[~, ~, ~, ~, ~, partsCoverArea, ~, ~, ~] = loadPartSelectionParameters(dataSetNumber);

% 1 + 1*1 + 2*2 + 3*3
totalNumBins = 1 + nbinsS1(1)*nbinsS1(2) + nbinsS2(1)*nbinsS2(2) + nbinsS3(1)*nbinsS3(2);


if numChannels == 1
    binSize = nNClusters{2};
elseif numChannels == 2
    binSize = nNClusters{2} + nNClusters{3};
elseif numChannels == 3    
    binSize = nNClusters{2} + nNClusters{3} + nNClusters{4};
end

%strRootE = 'D:/3D/Images for categorization/1view_1Scale/Elements';

strRootD = 'C:\Projects\Vladislav\Input data\Washington\Wash-rgbd-dataset_01_scale';

adder{2} = {'_layer2'};
adder{3} = {'_layer3'};
adder{4} = {'_layer4'};
adder{5} = {'_layer5'};
adder{6} = {'_layer6'};
adder{7} = {'_layer7'};
adder{8} = {'_layer8'};


is_subset = false;
subsetPercent = 1.0;

is_GPU_USED = true;


[~, list_mask, ~, lenF] = extractFileListWashington(strRootD, is_subset, subsetPercent);


for i = 1:numChannels 
    str = char(adder{i+1});
    str = [strRootD, str];
    [temp, model_id, category_id, lenF] = extractFileListWashingtonForClassification(str, is_subset, subsetPercent); % this folder with files of the second layer elements
    list_el{i} = temp;
end

addDescrs = 4;
descLength = addDescrs + binSize * totalNumBins;

tic

indd = randperm(lenF);

% input for the SVM
X = []; %zeros(len, descLength); 
Y = category_id(indd);                % class labels
M = model_id(indd);                % model number

parfor i = 1:lenF % for every image
    
    descriptor = zeros(1, descLength);

    I = [];
    addNNCl = 0;
    % address different viewing angles, scales
    for j = 1:numChannels 
         Itemp = imread(list_el{j}{indd(i)});
         Itemp = double(Itemp);
         Itemp = Itemp + addNNCl;
         Itemp(Itemp == addNNCl) = 0;
         I(:,:, j) = Itemp;
         addNNCl = addNNCl + nNClusters{j+1};
    end
    
    [r,c,ch] = size(I);  % elements
    mask = imread(list_mask{indd(i)});
    [rm, cm] = size(mask);
    maskArea = nnz(mask);
    
    if maskArea == 0 
        continue;
    end
    
    if rm ~= r || cm ~= c
        disp('Error');
    else
        [I, ~] = trimImage(I, mask);
    end
    
    [r,c,ch] = size(I);  % elements
    
    % compute the block size
    
    for sc = 1:3  % scales
        
        % 1*1*1 + 2*1*1 + 2*2*2 + 2*3*3
        if sc == 1
            nbins = nbinsS1;
            pos = binSize + 1; % after the first bin (overall)
        elseif sc == 2
            nbins = nbinsS2;
            pos = binSize * (1 + nbinsS1(1) * nbinsS1(2)) + 1;
        elseif sc == 3
            nbins = nbinsS3;
            pos = binSize * (1 + (nbinsS1(1) * nbinsS1(2) + nbinsS2(1) * nbinsS2(2))) + 1;
        end
    
        bx = floor(c/nbins(1)); % size of the bin in X direction
        by = floor(r/nbins(2)); % size of the bin in Y direction

        binOverall = zeros(1,binSize);

        for j = 1:nbins(1)
            for k = 1:nbins(2)

               binHist = zeros(1, binSize);  % current bin
               
               % compute boundaries of the bin
               x1 = (k-1)* bx + 1;
               x2 = k * bx;
               y1 = (j-1) * by + 1;
               y2 = j * by;
               
               % extract this bin from all images
               curI = I(y1:y2, x1:x2, :);
               curMask = mask(y1:y2, x1:x2, :);
               
               maskAreaLocal = nnz(curMask);
               if maskAreaLocal == 0
                   continue;
               end
               
               for ii = 1:numChannels

                    cover = partsCoverArea{ii + 1};

                    % create a structure element here
                    nm = [cover(2), cover(1)];
                    st = strel('rectangle', nm);


                    curIc = curI(:,:, ii);
                    
                    
                    while nnz(curIc) > 0
                        
                        maxEl = max(max(curIc));   
                        [rEls, cEls] = find(curIc == maxEl);
                        I2 = zeros(size(curIc));

                        indSS = sub2ind(size(I2), rEls, cEls); 
                        I2(indSS) = 1;
                        I2 = imdilate(I2, st);
                        curCover = nnz(I2);
                        binOverall(maxEl) = binOverall(maxEl) + curCover;
                        binHist(maxEl) = binHist(maxEl) + curCover;
                        
                        curIc(indSS) = 0;
                    end             
               end
               % normalize the histogram
               binHist = binHist/maskAreaLocal;

               % put binHist to the right position of the descriptor
               descriptor(pos:pos+binSize-1) = binHist;
               pos = pos + binSize;
            end
        end
    end
    
    % normalize the histogram binOverall
    binOverall = binOverall/maskArea;
    descriptor(1:binSize) = binOverall; 
    
    descriptor(end-3) = r;
    descriptor(end-2) = c;
    descriptor(end-1) = r/(r+c);
    descriptor(end) = maskArea/(r*c); 
    X = [X; descriptor];
    
    if mod(i,200) == 0
        i
    end
    
end
toc

save('xx4.mat', 'X', 'Y', 'M', '-v7.3');

end



