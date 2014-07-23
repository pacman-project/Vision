% this code is to create a histogram of second layer parts

% parameter with cross validation!
clear all;

nbinsS1 = [2,2];
nbinsS2 = [3,3];
totalNumBins = 1 + nbinsS1(1)*nbinsS1(2) + nbinsS2(1)*nbinsS2(2);
n2elements = 81;


%strRootE = 'D:/3D/Images for categorization/1view_1Scale/Elements';
strRootE = 'D:\3D\Images for categorization\1view_1Scale\Elements';
list_el = load_filelist(strRootE); % this folder with files of the second layer elements

histBin = zeros(1, n2elements);
len = length(list_el);

descLength = n2elements * totalNumBins;
descriptor = zeros(1, descLength);

% input for the SVM
X = zeros(len, descLength); 
Y = zeros(len, 1);                % class labels
M = zeros(len, 1);                % model number


for i = 1:len % for every image
    
    % address different viewing angles, scales
    I = imread(list_el{i});
    
    % find a class label first
    str = list_el{i};
    str = str(57:60);
    if str(2) == '_'
        ss = str(1);
        y = str2num(ss);
    elseif str(3) == '_'
        ss = str(1:2);
        y = str2num(ss);
    elseif str(4) == '_'
        ss = str(1:3);
        y = str2num(ss);
    end
    
    Y(i) = y;
    [r,c] = size(I); % second layer elements are in there
    % compute the block size
    
    for sc = 1:2
        
        if sc == 1
            nbins = nbinsS1;
            pos = 81 + 1;
        else
            nbins = nbinsS2;
            pos = 81 * 5 + 1;
        end
    
        bx = floor(c/nbins(1));
        by = floor(r/nbins(2));

        binOverall = zeros(1,81);

        for j =1:nbins(1)
            for k = 1:nbins(2)

               binHist = zeros(1,81); 
               % compute boundaries of the bin
               x1 = (k-1)* bx + 1;
               x2 = k * bx;
               y1 = (j-1) * by + 1;
               y2 = j * by;
               curI = I(y1:y2, x1:x2);

               els = curI(curI>0);
               ll = length(els);

               % now put all thes elements to the histogram
               for ii = 1:ll
                   curEl = els(ii);

                   binOverall(curEl) = binOverall(curEl) + 4;
                   binHist(curEl) = binHist(curEl) + 4;

                  % also increse score for the nighbour bins      
                   [ el_list, weights, numEl ] = extractNeighbours(curEl);
                   
                   for jj = 1:numEl
                       nEl = el_list(jj);
                       wei = weights(jj);
                       binOverall(nEl) = binOverall(nEl) + wei;
                       binHist(nEl) = binHist(nEl) + wei;
                   end                            
               end

               % normalize the histogram
               maxH = max(binHist);
               binHist = binHist/maxH;

               % put binHist to the right position of the descriptor
               descriptor(pos:pos+80) = binHist;
               pos = pos + 81;

            end
        end
    end
    
    % normalize the histogram binOverall
    maxH = max(binOverall);
    binOverall = binOverall/maxH;
    descriptor(1:81) = binOverall;   
    X(i,:) = descriptor;
    
    if mod(i,20) == 0
        i
    end
    
end

M = Y;
Y = floor((Y - 1) / 20) + 1;
save('xx.mat', 'X', 'Y', 'M');





