% this code is to create a histogram of compositional parts
function [] = histogramOfParts_best(n2Clusters, n3Clusters, n4Clusters)

nbinsS1 = [1,1];
nbinsS2 = [2,2];
nbinsS3 = [3,3];

are_neighboursInvolved = false;

% 1 + 1*1 + 2*2 + 3*3
totalNumBins = 1 + nbinsS1(1)*nbinsS1(2) + nbinsS2(1)*nbinsS2(2) + nbinsS3(1)*nbinsS3(2);

numChannels = 3;

if numChannels == 1
    binSize = n2Clusters;
elseif numChannels == 2
    binSize = n2Clusters + n3Clusters;
elseif numChannels == 3    
    binSize = n2Clusters + n3Clusters + n4Clusters;
end

%strRootE = 'D:/3D/Images for categorization/1view_1Scale/Elements';
strRootE = 'D:\3D\Input Data\Images for categorization\8view_1Scale\ElementsOpt3';
strRootD = 'D:\3D\Input Data\Images for categorization\8view_1Scale\images';

list_el = load_filelist(strRootE); % this folder with files of the second layer elements
list_D = load_filelist(strRootD); % this is folder with depths
len = length(list_el);

descLength = binSize * totalNumBins;

% input for the SVM
X = [];%zeros(len, descLength); 
Y = zeros(len, 1);                % class labels
M = zeros(len, 1);                % model number

weightsAll = zeros(n2Clusters,8);
el_listAll = zeros(n2Clusters,8);
numElsAll = zeros(n2Clusters);
 
if are_neighboursInvolved  % precompute the table for speed up reasons
    for i = 1:n2Clusters
        [ el_list, weights, numEl ] = extractNeighbours(i);
        el_listAll(i, :) = el_list;
        weightsAll(i, :) = weights;
        numElsAll(i) = numEl;
    end
end

parfor i = 1:len % for every image
    
    descriptor = zeros(1, descLength);
    % address different viewing angles, scales
    I = imread(list_el{i});
    
    I_D = imread(list_D{i});
    [r,c] = size(I_D);    % elements
    I_D = double(I_D);
    
    H = I_D(:);
    H = H(H>0);
    
    % trim image I_D and chech whether it's size is the same as size of I
    mask = zeros(r,c);
    mask(I_D > 0) = 1;
    [I_D, ~] = trimImage(I_D, mask);
    [r,c] = size(I_D);
    [r1, c1, ~] = size(I); 
    
    if r ~= r1 || c ~=c1
        disp('ERROR');  % should not happen
    end
    
    % find a class label first
    str = list_el{i};
    index_ = strfind(str,'/');
    new_str= str(index_+1:end);
    index = strfind(new_str,'_');
    no = new_str(index(1)+1:index(2)-1);
    y = str2num(no);
    
%     str = str(57:60);
%     if str(2) == '_'
%         ss = str(1);
%         y = str2num(ss);
%     elseif str(3) == '_'
%         ss = str(1:2);
%         y = str2num(ss);
%     elseif str(4) == '_'
%         ss = str(1:3);
%         y = str2num(ss);
%     end
%     
    Y(i) = y;
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
               curI_D = I_D(y1:y2, x1:x2);
               
               for ii = 1:numChannels

                   curIc = curI(:,:, ii)
                   % extract all non-zero elements
                   [rEls,cEls] = find(curIc > 0);
                   ll = length(rEls);

                   % now put all these elements to the histogram
                   for ii = 1:ll
                       curEl = curIc(rEls(ii), cEls(ii));    
                       binOverall(curEl) = binOverall(curEl) + 4;  % feed the element to the overall bin
                       binHist(curEl) = binHist(curEl) + 4; 

                      % also increase score for the nighbour bins 
                      if are_neighboursInvolved           
                           numEl = numElsAll(curEl);

                           for jj = 1:numEl
                               nEl = el_listAll(curEl, jj);
                               wei = weightsAll(curEl, jj);
                               binOverall(nEl) = binOverall(nEl) + wei;
                               binHist(nEl) = binHist(nEl) + wei;
                           end   
                      end
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
    X = [X; descriptor];
    
    if mod(i,20) == 0
        i
    end
    
end

M = Y;
Y = floor((Y - 1) / 20) + 1;
save('xx.mat', 'X', 'Y', 'M');

end





