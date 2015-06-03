% this function computes coverage of all layers and writes the result to
% the folder

function [] = performInference2(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, nClusters, ...
                            infArray, outRoot, combs, largestLine, wCoverage, wOverlap, is_inhibition, ...
                            is_downsampling, dowsample_rate, dataSetNumber, depthPath, cluster1Centres, cluster1Bounds, ...
                            is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, isReconstructionError, ...
                            dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward, is_overwrite)
            

% -------------for the preliminary processing function--------------------
isY = true;
isX = true;
isX_FB = true;
isTrim = false;

discrThresh1 = 10;
discrThresh2 = 3;

lenDPW = length(depthPath);

if dataSetNumber ~= 2
    list_mask = zeros(1, lenF);
end

indsSS = randperm(lenF);

parfor i = 1:lenF
    
    curStr = list_depth{indsSS(i)};
%       ll = strfind(curStr, '/');
%       fileName = curStr(ll:end);
        

    fileName = curStr(lenDPW+1:end);
    outFile = [outRoot, fileName];
    
    ll = strfind(outFile, '/');
    lll = strfind(outFile, '\');
    ll = [ll, lll];
    ll = max(ll); % last position
    folderName = outFile(1:ll);
    b = exist(folderName,'dir');

    if b == 0
        mkdir(folderName);
    end
    % chech if file exists
    b = exist(outFile, 'file');

   if ~b || is_overwrite
    
    I = imread(list_depth{indsSS(i)}); % list should contain the full path
    I = I(:,:,1);
    I = double(I);
    mask = [];
    
    if dataSetNumber == 1 || dataSetNumber == 3 || dataSetNumber == 4    % Aim@Shape dataset
        
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discSize, isX, isY, isX_FB, ...
                            isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                            dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward);
                        
    
    elseif dataSetNumber == 2  % Washington data set
        
        mask = imread(list_mask{indsSS(i)});

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
              
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, isX, isY, isX_FB,...
                    isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                    dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward);
                
    end
    
    if ~is_successfull
        continue;
    end 
    
    marks1 = zeros(r,c);       % here we collect marks
    marks1_alt = zeros(r,c);   % this are second choice marks (with larger error)
    errors1 = zeros(r,c);      % here we collect errors for each mark
    errors1_alt = zeros(r,c);  % here we collect errors for each mark
    
    
%% RECOGNITION OF THE FIRST LAYER
    
    % perform line discretization in X direction (TAKE ODD LINES ONLY)
    
    for k = 1:2:r % for every other line we perform discretization
        ind = find(mask(k,:));
        if length(ind) < discrThresh1
            marks1(k,:) = zeros(1, c); % line of zeros in this case
        else
            % here is discretization itself
            if mod(ind(1),2) == 1 % we always start from numbers such that: number mod 2 = 0
                offset = 1;
            elseif mod(ind(1),2) == 0
                offset = 0;
            end
            
            
            for iii = 1:3
                
                fx = Ix(k,ind(1)+offset : ind(end), iii);  % look to central difference
                strLen = length(fx);
                
                if iii == 1
                    outLine = zeros(1,strLen);
                    outErrors = zeros(1,strLen);
                    altClusters = zeros(1,strLen);
                    altErrors = zeros(1,strLen);
                end
                
                inds = 1:2:strLen-1;
                fxQuant = fx(inds);
                strLen = length(fxQuant);
                
                if iii == 1
                    [nearestClusters,  alternativeClusters] = discretizeLineVectorized(fxQuant, nClusters, cluster1Centres, cluster1Bounds);
%                 elseif iii == 2  % check if forward and backward differences are the same
%                     [nearestClusters2,  alternativeClusters2] = discretizeLineVectorized(fxQuant, nClusters, cluster1Centres, cluster1Bounds);
%                     nearestClusters(nearestClusters ~= nearestClusters2) = 0; % non-planar area!
                end
 
            end
            
            
            % [errors, alternativeErrors] = computeLineReconstructionErrors(nearestClusters, alternativeClusters, cluster1Centres, fx, inds, thresh);
            
            if is_inhibition{1}
                [output, curErrs] = lineInhibition(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine);
                outLine(inds) = output;
                outErrors(inds) = curErrs;
                % we do not perform inhibition in alternative hypothesis space
            else
                outLine(inds) = nearestClusters;
            end
            altClusters(inds) = alternativeClusters;
%             altErrors(inds) = alternativeErrors;
            
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
    
    if infArray(1)  % write the outcome somewhere
%            imtool(marks1, [0,nClusters]);
%           imtool(I, [min(min(I)), max(max(I))]);
%         imtool(errors1);
%         imtool(marks1_alt, [0,9]);
%         imtool(errors1_alt);
%         TODO write the results somewhere
    end
    
    if isReconstructionError{1}
        disp('not implemented yet');
    end
    

    a = 2;
    
    
%% RECOGNITION OF THE SECOND LAYER 
    
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
    
    % perform line discretization in Y direction (TAKE columns such that number mod 2 = 0)
    
    for k = 4:2:c-2 % for every other column we perform discretization
        
        ind = find(marks1(:,k));
        if length(ind) < discrThresh2
            marks2(:,k) = zeros(r, 1); % line of zeros in this case
        else     
            
            % extract slices in y-direction
            yMarks1 = marks1(ind(1) : ind(end), k-2: k+2);
            yErrors1 = errors1(ind(1) : ind(end), k-2: k+2);
            yMarks1_alt = marks1_alt(ind(1) : ind(end), k-2: k+2);
            yErrors1_alt = errors1_alt(ind(1) : ind(end), k-2: k+2);
            
            line_Iy1 = Iy(ind(1) : ind(end), k, 1);
            line_Iy2 = Iy(ind(1) : ind(end), k, 2);
            line_Iy3 = Iy(ind(1) : ind(end), k, 3);
            
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
                dY = line_Iy1(indsL(jj));
                clusterY = define1Cluster(dY, cluster1Bounds, nClusters); 
                
                dY2 = line_Iy2(indsL(jj));
                clusterY2 = define1Cluster(dY2, cluster1Bounds, nClusters);
                dY3 = line_Iy3(indsL(jj));
                clusterY3 = define1Cluster(dY3, cluster1Bounds, nClusters);
                
                if clusterY == 0 || clusterY~=clusterY2 || clusterY~=clusterY3
                    nearestClusters(indsL(jj)) = 0;
                else
                    errorY = 0.1;                                                           % TODO SOMETHING BETTER
                    newInd = compute2elementIndex(clusterX, clusterY, nClusters);
                    nearestClusters(indsL(jj)) = newInd;
                    errors(indsL(jj)) = (errors(indsL(jj)) + errorY)/2;  % controversal
                end
            end
            
            if is_inhibition{2} % even numbers
                inds = 3:2:strLen;
                nearestClusters = nearestClusters(inds);
                errors = errors(inds);
                [output, curErrs] = lineInhibition(nearestClusters', errors', wCoverage, wOverlap, combs, largestLine);    
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
        
    if infArray(2) % write the results somewhere
        
%         imtool(I, [min(min(I)), max(max(I))]);
%         imtool(Ix, [min(min(Ix)), max(max(Ix))]);
%         imtool(Iy, [min(min(Iy)), max(max(Iy))]);
%         imtool(marks2, [0, nClusters^2]);    % GIVES NICE PICTURES

        marks2 = uint8(marks2);

        imwrite(marks2, outFile, 'png');
        
    end
    
    if isReconstructionError{2}
        disp('not implemented');
    end
    
    

    if mod(i,10) == 0
        i
    end
   end
        
end












