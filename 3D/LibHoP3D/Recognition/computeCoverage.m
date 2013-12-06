% this function computes coverage of all layers and writes the result to
% the folder

function [] = computeCoverage(list_depth, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, cluster1Centres, cluster1Lengths, thresh, ...
                areLayersRequired, outRoot, combs, largestLine, wCoverage, wOverlap)
 
% -------------for the preliminary processing function--------------------
isErrosion = false;
discSize = 0;
isY = true;
isTrim = true;

discrThresh1 = 10;
discrThresh2 = 2;

for i = 1:1  % To use later for models
    
    I = imread(list_depth{i}); % list should contain the full path
    %imtool(I);
    [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, isErrosion, discSize, isY, isTrim, dxKernel, sigmaKernelSize, sigma);

    if ~is_successfull
        continue;
    end 
    marks1 = zeros(r,c); % here we collect marks
    errors1 = zeros(r,c); % here we collect errors for each mark
    
    % perform line discretization in X direction (TAKE ODD LINES ONLY)
    for k = 1:1:r % for every line we perform a discretization
        ind = find(mask(k,:));
        if length(ind) < discrThresh1
            marks1(k,:) = zeros(1, c); % line of zeros in this case
        else
            % here is discretization itself
            if mod(ind(1),2) ~= 0 % we always start from even numbers
                offset = 1;
            else
                offset = 0;
            end
                
            fx = Ix(k,ind(1)+offset : ind(end));
            strLen = length(fx);
            outLine = zeros(1,strLen);
            outErrors = zeros(1,strLen);
            inds = 1:2:strLen-1;
            fx = fx(inds);
            strLen = length(fx);
            [nearestClusters, errors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, cluster1Lengths, thresh); % replaces filter responces with nearest clusters
            [output, curErrs] = lineDiscretizationOptLayer1(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine);     
            outLine(inds) = output;
            outErrors(inds) = curErrs;
            marks1(k,ind(1)+offset:ind(end)) = outLine;
            errors1(k,ind(1)+offset:ind(end)) = outErrors;
        end
    end

    marks1 = marks1 .* mask;
    errors1 = errors1 .* mask;
    imtool(marks1, [0,9]);
%     imtool(errors1);
    
%   I(marks1>0) = 65000;
%   I = uint16(I);
%   imtool(I,[]);


    % here we consider the SECOND LAYER COVERAGE (THIS IS YET A PROTOTYPED VERSION)
    
    % add zero boundaries to the marks1 end errors1 images
    offset2 = 2;
    
    marks1 = addZeroBoundaries(marks1, offset2);
    errors1 = addZeroBoundaries(errors1, offset2);
    
    [r,c] = size(marks1);
    marks2 = zeros(r,c);
    marksTrial = zeros(r,c);
    errors2 = zeros(r,c);
    
    % perform line discretization in Y direction (TAKE EVEN COLUMNS ONLY!)
    for k = 4:2:c % for every column we perform discretization
        ind = find(marks1(:,k));
        if length(ind) < discrThresh2
            marks2(:,k) = zeros(r, 1); % line of zeros in this case
        else         
            % extract slices
            sliceMarks1 = marks1(ind(1) : ind(end), k-2:k+2);
            sliceErrors1 = errors1(ind(1) : ind(end), k-2:k+2);
            
            [strLen,~] = size(sliceMarks1);
            
            outLine = zeros(strLen, 1);
            outErrors = zeros(strLen, 1);
            
            [nearestClusters, errors] = discretizeSlice(sliceMarks1, sliceErrors1, strLen); % prepare data for line discretization
            marksTrial(ind(1):ind(end), k) = nearestClusters';
            inds = 2:2:strLen;
            nearestClusters = nearestClusters(inds);
            errors = errors(inds);
            [output, curErrs] = lineDiscretizationOptLayer1(nearestClusters', errors', wCoverage, wOverlap, combs, largestLine);    
            outLine(inds) = output';
            outErrors(inds) = curErrs';
            marks2(ind(1):ind(end), k) = outLine;
            errors2(ind(1):ind(end), k) = outErrors;
        end
    end
    
    % return to the original size
    marksTrial = marksTrial(offset2+1 : r-offset2, offset2+1 : c-offset2);
    marks2 = marks2(offset2+1 : r-offset2, offset2+1 : c-offset2);
    errors2 = errors2(offset2+1 : r-offset2, offset2+1 : c-offset2);
    
    imtool(marks2, [0, 9]);  %   GIVES NICE PICTURES
    imtool(marksTrial, [0, 9]);  %   GIVES NICE PICTURES
    
    % make them second layer elements
    [rows, cols] = find(marks2 > 0);
    ll = length(rows);
    for jj = 1:ll
        dY = Iy(rows(jj), cols(jj));
        clusterY = define1Cluster(dY, nClusters, cluster1Lengths, thresh);
        newInd = compute2elementIndex(marks2(rows(jj), cols(jj)), clusterY, nClusters);
        marks2(rows(jj), cols(jj)) = newInd;
    end
   
    
%    imtool(marks2, [0, 81]);  %   GIVES NICE PICTURES       

% 
%     % now all marks are in range 1..81      
%     % write them to the file        
%     marks = uint8(marks);
%     curStr = list_depth{i};
%     filename = curStr(sl+2:end);
%     outFile = [strOutput, filename];
%     imwrite(marks, outFile, 'png');
% 
% %         IG = uint16(IG);
% %         outFile = [strOutput, strOutTI, str, '.png'];
% %         imwrite(IG, outFile, 'png');
% 
    if mod(i,10) == 0
        i
    end
        
end



% % This is the old code fragment for the second layer detection
%     for jj = offset2+1 : r-offset2 
%         for kk = offset2+1 : c-offset2
%             central = marks1(jj,kk);
%             if central ~= 0
%                 window = marks1(jj-2:jj+2, kk-2:kk+2);
%                 kol = length(find(window == central));
%                 if central ~= 1
%                     kol = kol + 0.5 * length(find(window == central - 1)); % we can do something smart with errors here
%                 end
%                 if central ~= nClusters
%                     kol = kol + 0.5 * length(find(window == central + 1));
%                 end
%                 if kol >= 4
%                     %  take a Y derivative in central point and define a
%                     %  cluster
%                     dY = Iy(jj,kk); 
%                     clusterY = define1Cluster(dY, nClusters, cluster1Lengths, thresh);                      
%                     marks2(jj,kk) = compute2elementIndex(central, clusterY, nClusters);
%                     window (window == central) = 0;
%                     marks1(jj-2:jj+2, kk-2:kk+2) = window;
% 
% %                     window = window(2:4,2:4);
% %                     window(window == central) = 0;
% %                     marks1(jj-1:jj+1, kk-1:kk+1) = window;
%                     
%                 end
%             end
%         end
%     end












