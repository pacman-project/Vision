% this function computes coverage of all layers and writes the result to
% the folder

function [] = computeCoverage(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, nClusters, ...
                            vocabulary1Layer, areLayersRequired, outRoot, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                            is_downsampling, dowsample_rate, dataSetNumber, depthPath, cluster1Centres, cluster1Lengths, thresh)
            

% -------------for the preliminary processing function--------------------
isY = true;
isX = true;
isTrim = true;

discrThresh1 = 10;
discrThresh2 = 3;

lenDPW = length(depthPath);

is_first_layer = areLayersRequired(1);
is_second_layer = areLayersRequired(2);
is_inhibition_1_layer = isInhibitionRequired(1);
is_inhibition_2_layer = isInhibitionRequired(2); 

parfor i = 1:lenF  % To use later for models
    
    I = imread(list_depth{i}); % list should contain the full path
    I = I(:,:,1);
    I = double(I);
    mask = [];
    
    if dataSetNumber == 1      % Aim@Shape dataset
        
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discSize, isX, isY, ...
                            isTrim, dxKernel, sigmaKernelSize, sigma);
    
    elseif dataSetNumber == 2  % Washington data set
        
        mask = imread(list_mask{i});

        if is_downsampling

            [r,c] = size(I);
            maxSize = max(r,c);
            downRate = dowsample_rate;

            if maxSize > 100                   % controversal rescaling method
                downRate = 250 /  maxSize;
            end
            if downRate < 1
                downRate = 1;
            end

            I = imresize(I, downRate);
            mask = imresize(mask, downRate);
        end
        
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, isX, isY,...
                    isTrim, dxKernel, sigmaKernelSize, sigma);        
    end
    
    if ~is_successfull
        continue;
    end 
    
    marks1 = zeros(r,c);       % here we collect marks
    marks1_alt = zeros(r,c);   % this are second choice marks (with larger error)
    errors1 = zeros(r,c);      % here we collect errors for each mark
    errors1_alt = zeros(r,c);  % here we collect errors for each mark
    
    
    % RECOGNITION OF THE FIRST LAYER --------------------------------------
    
    % perform line discretization in X direction (TAKE ODD LINES ONLY)
    
    for k = 1:2:r % for every other line we perform discretization
        ind = find(mask(k,:));
        if length(ind) < discrThresh1
            marks1(k,:) = zeros(1, c); % line of zeros in this case
        else
            % here is discretization itself
            if mod(ind(1),3) == 2 % we always start from numbers such that: number mod 3 = 0
                offset = 1;
            elseif mod(ind(1),3) == 1
                offset = 2;  
            elseif mod(ind(1),3) == 0
                offset = 0;
            end
                
            fx = Ix(k,ind(1)+offset : ind(end));
            strLen = length(fx);
            outLine = zeros(1,strLen);
            outErrors = zeros(1,strLen);
            altClusters = zeros(1,strLen);
            altErrors = zeros(1,strLen);
            inds = 1:3:strLen-1;
            fx = fx(inds);
            strLen = length(fx);
            

            [nearestClusters, errors, alternativeClusters, alternativeErrors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, ...
                cluster1Lengths, thresh); % replaces filter responces with nearest clusters
            
            if is_inhibition_1_layer
                [output, curErrs] = lineDiscretizationOptLayer1(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine);
                outLine(inds) = output;
                outErrors(inds) = curErrs;
                % we do not perform inhibition in alternative hypothesis space
            else
                outLine(inds) = nearestClusters;
                outErrors(inds) = errors;
            end
            altClusters(inds) = alternativeClusters;
            altErrors(inds) = alternativeErrors;
            
            marks1(k,ind(1)+offset:ind(end)) = outLine;
            errors1(k,ind(1)+offset:ind(end)) = outErrors;
            marks1_alt(k,ind(1)+offset:ind(end)) = altClusters;
            errors1_alt(k,ind(1)+offset:ind(end)) = altErrors;
        end
    end

    marks1 = marks1 .* mask;
    errors1 = errors1 .* mask;
    marks1_alt = marks1_alt .* mask;
    errors1_alt = errors1_alt.*mask;
    
    if is_first_layer  % write the outcome somewhere
%          imtool(marks1, [0,9]);
%         imtool(errors1);
%         imtool(marks1_alt, [0,9]);
%         imtool(errors1_alt);
%         TODO write the results somewhere
    end

    a = 2;
    
    
    % RECOGNITION OF THE SECOND LAYER --------------------------------------
    
    % add zero boundaries to the marks1 end errors1 images
%     offset2 = 2;
%     
%     marks1 = addZeroBoundaries(marks1, offset2);
%     errors1 = addZeroBoundaries(errors1, offset2);
%     marks_alt = addZeroBoundaries(marks1_alt, offset2);
%     errors_alt = addZeroBoundaries(errors1_alt, offset2);
    
    [r,c] = size(marks1);
    marks2 = zeros(r,c);
    errors2 = zeros(r,c);
    
    % perform line discretization in Y direction (TAKE columns such that number mod 3 = 0)
    for k = 3:3:c % for every column we perform discretization
        
        ind = find(marks1(:,k));
        if length(ind) < discrThresh2
            marks2(:,k) = zeros(r, 1); % line of zeros in this case
        else         
            % extract lines in y-direction
            yMarks1 = marks1(ind(1) : ind(end), k);
            yErrors1 = errors1(ind(1) : ind(end), k);
            yMarks1_alt = marks1_alt(ind(1) : ind(end), k);
            yErrors1_alt = errors1_alt(ind(1) : ind(end), k);
            line_Iy = Iy(ind(1) : ind(end), k);
            
            if mod(ind(1), 2) == 0
                offset = 1;   
            elseif mod(ind(1), 2) == 1  % should always be the case!
                offset = 0;  
            end
            
            strLen = size(yMarks1, 1);
            outLine = zeros(strLen, 1);
            outErrors = zeros(strLen, 1);
            [nearestClusters, errors] = discretizeY(yMarks1, yErrors1, yMarks1_alt, yErrors1_alt, strLen, offset); % prepare data for line discretization
            
            % make them second layer elements
            indsL = find(nearestClusters > 0);
            for jj = 1:length(indsL)
                clusterX = nearestClusters(indsL(jj));
                dY = line_Iy(indsL(jj));
                clusterY = define1Cluster(dY, nClusters, cluster1Lengths, thresh); 
                errorY = abs(  (dY - cluster1Centres(clusterY)) / cluster1Lengths(clusterY)  ); % compute an error
                if errorY > 1
                    errorY = 1;
                end
                newInd = compute2elementIndex(clusterX, clusterY, nClusters);
                nearestClusters(indsL(jj)) = newInd;
                errors(indsL(jj)) = (errors(indsL(jj)) + errorY)/2;  % controversal
            end
            
            if is_inhibition_2_layer  % even numbers
                inds = 2:2:strLen;
                nearestClusters = nearestClusters(inds);
                errors = errors(inds);
                [output, curErrs] = lineDiscretizationOptLayer1(nearestClusters', errors', wCoverage, wOverlap, combs, largestLine);    
                outLine(inds) = output';
                outErrors(inds) = curErrs';
                marks2(ind(1):ind(end), k) = outLine;
                errors2(ind(1):ind(end), k) = outErrors;
            else
                marks2(ind(1):ind(end), k) = nearestClusters';
            end
        end
    end
    % return to the original size
    
%     marksTrial = marksTrial(offset2+1 : r-offset2, offset2+1 : c-offset2);
%     marks2 = marks2(offset2+1 : r-offset2, offset2+1 : c-offset2);
%     errors2 = errors2(offset2+1 : r-offset2, offset2+1 : c-offset2);
        
    if is_second_layer % write the results somewhere
        
%          imtool(I, [min(min(I)), max(max(I))]);
%         imtool(Ix, [min(min(Ix)), max(max(Ix))]);
%         imtool(Iy, [min(min(Iy)), max(max(Iy))]);
%          imtool(marks2, [0, 81]);    % GIVES NICE PICTURES

        marks2 = uint8(marks2);
        curStr = list_depth{i};
%       ll = strfind(curStr, '/');
%       fileName = curStr(ll:end);

        fileName = curStr(lenDPW+1:end);
        outFile = [outRoot, fileName];
        
        ll = strfind(outFile, '/');
        ll = ll(end); % last position
        folderName = outFile(1:ll);
        b = exist(folderName,'dir');
        
        if b == 0
            mkdir(folderName);
        end
        imwrite(marks2, outFile, 'png');
        
%         IG = uint16(IG);
%         outFile = [strOutput, strOutTI, str, '.png'];
%         imwrite(IG, outFile, 'png');
    end
    
    
   % RECOGNITION OF THE THIRD LAYER --------------------------------------

    if mod(i,100) == 0
        i
    end

        
end












