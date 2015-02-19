% this code is to create a MULTISCALE histogram of compositional parts

function [] = histogramOfParts_best_Washington(nNClusters, numChannels)

% pyramid_levels = [1.0, 0,66];
% pyramid_num_levels = length(pyramid_levels);

nbinsS1 = [1,1];
nbinsS2 = [2,2];
nbinsS3 = [3,3];

are_neighboursInvolved = false;

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

strRootD = 'D:\Input Data\Washington\Wash-rgbd-dataset_02_scale';

adder{2} = {'_layer2'};
adder{3} = {'_layer3'};
adder{4} = {'_layer4'};
adder{5} = {'_layer5'};
adder{6} = {'_layer6'};
adder{7} = {'_layer7'};
adder{8} = {'_layer8'};


is_subset = false;
subsetPercent = 1.0;


[~, list_mask, ~, lenF] = extractFileListWashington(false, strRootD, '', is_subset, subsetPercent);


for i = 1:numChannels 
    str = char(adder{i+1});
    str = [strRootD, str];
    [temp, model_id, category_id, lenF] = extractFileListWashingtonForClassification(str, is_subset, subsetPercent); % this folder with files of the second layer elements
    list_el{i} = temp;
end

addDescrs = 4;
descLength = addDescrs + binSize * totalNumBins;

% input for the SVM
X = [];%zeros(len, descLength); 
Y = category_id;                % class labels
M = model_id;                % model number

 
% if are_neighboursInvolved  % precompute the table for speed up reasons
%     for i = 1:nNClusters{2}
%         [ el_list, weights, numEl ] = extractNeighbours(i);
%         el_listAll(i, :) = el_list;
%         weightsAll(i, :) = weights;
%         numElsAll(i) = numEl;
%     end
% end

parfor i = 1:lenF % for every image
    
    descriptor = zeros(1, descLength);

    I = [];
    addNNCl = 0;
    % address different viewing angles, scales
    for j = 1:numChannels 
         Itemp = imread(list_el{j}{i});
         Itemp = double(Itemp);
         Itemp = Itemp + addNNCl;
         Itemp(Itemp == addNNCl) = 0;
         I(:,:, j) = Itemp;
         addNNCl = addNNCl + nNClusters{j+1};
    end
    
    [r,c,ch] = size(I);  % elements
    mask = imread(list_mask{i});
    [rm, cm] = size(mask);
    maskArea = length(mask(mask == 1));
    
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
               
               % extract this bin from both images
               curI = I(y1:y2, x1:x2, :);
               
               for ii = 1:numChannels

                   curIc = curI(:,:, ii);
                   % extract all non-zero elements
                   [rEls,cEls] = find(curIc > 0);
                   ll = length(rEls);

                   % now put all these elements to the histogram
                   for iii = 1:ll
                       curEl = curIc(rEls(iii), cEls(iii));    
                       binOverall(curEl) = binOverall(curEl) + 4;  % feed the element to the overall bin
                       binHist(curEl) = binHist(curEl) + 4; 

                      % also increase score for the nighbour bins 
%                       if are_neighboursInvolved           
%                            numEl = numElsAll(curEl);
% 
%                            for jj = 1:numEl
%                                nEl = el_listAll(curEl, jj);
%                                wei = weightsAll(curEl, jj);
%                                binOverall(nEl) = binOverall(nEl) + wei;
%                                binHist(nEl) = binHist(nEl) + wei;
%                            end   
%                       end
                   end
                   
               end

               % normalize the histogram
               maxH = max(binHist);
               if maxH ~= 0
                    binHist = binHist/maxH;
               end
               % put binHist to the right position of the descriptor
               descriptor(pos:pos+binSize-1) = binHist;
               pos = pos + binSize;
            end
        end
    end
    
    % normalize the histogram binOverall
    maxH = max(binOverall);
    binOverall = binOverall/maxH;
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

save('xx4.mat', 'X', 'Y', 'M', '-v7.3');

end




% % this code is to create a MULTISCALE histogram of compositional parts
% 
% function [] = histogramOfParts_best_Washington(numChannels, nNClusters)
% 
% % nNClusters{2} = 81;
% % nNClusters{3} = 500;
% % nNClusters{4} = 700;
% 
% disp('Building a histogram of parts...');
% 
% nbinsS1 = [1,1];
% nbinsS2 = [2,2];
% nbinsS3 = [3,3];
% 
% % are_neighboursInvolved = false;
% 
% % 1 + 1*1 + 2*2 + 3*3
% totalNumBins = 1 + nbinsS1(1)*nbinsS1(2) + nbinsS2(1)*nbinsS2(2) + nbinsS3(1)*nbinsS3(2);
% 
% 
% 
% if numChannels == 1
%     binSize = nNClusters{2};
% elseif numChannels == 2
%     binSize = nNClusters{2} + nNClusters{3};
% elseif numChannels == 3    
%     binSize = nNClusters{2} + nNClusters{3} + nNClusters{4};
% end
% 
% adder{2} = {'_layer2'};
% adder{3} = {'_layer3'};
% adder{4} = {'_layer4'};
% adder{5} = {'_layer5'};
% adder{6} = {'_layer6'};
% adder{7} = {'_layer7'};
% adder{8} = {'_layer8'};
% 
% strRootD = 'D:\Input Data\Washington\Old Files\Wash-rgbd-dataset_007';
% is_subset = false;
% subsetPercent = 1.0;
% 
% % [~, list_mask, ~, lenF] = extractFileListWashington(false, strRootD, [], is_subset, subsetPercent); % list of masks
% 
% for i = 1:numChannels 
%     str = char(adder{i+1});
%     str = [strRootD, str];
%     [temp, model_id, category_id, lenF] = extractFileListWashingtonForClassification(str, is_subset, subsetPercent); % this folder with files of the second layer elements
%     list_el{i} = temp;
% end
% 
% genDescs = 0;
% 
% descLength = genDescs + binSize * totalNumBins;
% 
% % input for the SVM
% X = [];                % descriptors 
% Y = category_id;
% M = model_id;
% 
%  
% % if are_neighboursInvolved  % precompute the table for speed up reasons
% %     for i = 1:nNClusters{2}
% %         [ el_list, weights, numEl ] = extractNeighbours(i);
% %         el_listAll(i, :) = el_list;
% %         weightsAll(i, :) = weights;
% %         numElsAll(i) = numEl;
% %     end
% % end
% 
% 
% parfor i = 1:lenF % for every image
%     
%     descriptor = zeros(1, descLength);
%     
%     I = [];
%     
%     addNNCl = 0;
% 
%     for j = 1:numChannels  % combine all channels to one image
%          Itemp = imread(list_el{j}{i});
%          if j == 1
%              [r,c] = size(Itemp);
%              I = zeros(r,c,numChannels);
%          end
%          Itemp = double(Itemp);
%          Itemp = Itemp + addNNCl;
%          Itemp(Itemp == addNNCl) = 0;
%          I(:,:, j) = Itemp;
%          addNNCl = addNNCl + nNClusters{j+1};
%     end
%     
% %     mask = imread(list_mask{i});
% %     [rm, cm] = size(mask);
% %     [r,c,ch] = size(I);  % elements
% %     
% %     if (rm ~= r) || (cm~=c)
% %         disp('Error');
% %         return;
% %     else
% %         % trim image according to it's mask
% %         [I, ~] = trimImage(I, mask);
% %     end
%     
%     [r,c,ch] = size(I);  % elements
%     
%     % compute the block size
%     
%     for sc = 1:3  % scales
%         
%         % 1*1*1 + 2*1*1 + 2*2*2 + 2*3*3
%         if sc == 1
%             nbins = nbinsS1;
%             pos = binSize + 1; % after the first bin (overall)
%         elseif sc == 2
%             nbins = nbinsS2;
%             pos = binSize * (1 + nbinsS1(1) * nbinsS1(2)) + 1;
%         elseif sc == 3
%             nbins = nbinsS3;
%             pos = binSize * (1 + (nbinsS1(1) * nbinsS1(2) + nbinsS2(1) * nbinsS2(2))) + 1;
%         end
%     
%         bx = floor(c/nbins(1)); % size of the bin in X direction
%         by = floor(r/nbins(2)); % size of the bin in Y direction
% 
%         binOverall = zeros(1,binSize);
% 
%         for j = 1:nbins(1)
%             for k = 1:nbins(2)
% 
%                binHist = zeros(1, binSize);  % current bin
%                
%                % compute boundaries of the bin
%                x1 = (k-1)* bx + 1;
%                x2 = k * bx;
%                y1 = (j-1) * by + 1;
%                y2 = j * by;
%                
%                % extract this bin from both images
%                curI = I(y1:y2, x1:x2, :);
%                
%                for ii = 1:numChannels
% 
%                    curIc = curI(:,:, ii);
%                    % extract all non-zero elements
%                    [rEls,cEls] = find(curIc > 0);
%                    ll = length(rEls);
% 
%                    % now put all these elements to the histogram
%                    for ii = 1:ll
%                        curEl = curIc(rEls(ii), cEls(ii));    
%                        binOverall(curEl) = binOverall(curEl) + 4;  % feed the element to the overall bin
%                        binHist(curEl) = binHist(curEl) + 4; 
% 
%                       % also increase score for the nighbour bins 
% %                       if are_neighboursInvolved           
% %                            numEl = numElsAll(curEl);
% % 
% %                            for jj = 1:numEl
% %                                nEl = el_listAll(curEl, jj);
% %                                wei = weightsAll(curEl, jj);
% %                                binOverall(nEl) = binOverall(nEl) + wei;
% %                                binHist(nEl) = binHist(nEl) + wei;
% %                            end   
% %                       end
%                    end
%                    
%                end
% 
%                % normalize the histogram
%                maxH = max(binHist);
%                if maxH ~= 0
%                     binHist = binHist/maxH;
%                end
%                % put binHist to the right position of the descriptor
%                descriptor(pos:pos+binSize-1) = binHist;
%                pos = pos + binSize;
%             end
%         end
%     end
%     
%     % normalize the histogram binOverall
%     maxH = max(binOverall);
%     binOverall = binOverall/maxH;
%     descriptor(1:binSize) = binOverall;
%     
% %     descriptor(end - 1) = 0;  % no additional features yet
% %     descriptor(end) = 0;
%     
%     X = [X; descriptor];
%     
%     if mod(i,20) == 0
%         i
%     end
%     
% end
% 
% save('xx4.mat', 'X', 'Y', 'M');
% 
% end
% 




