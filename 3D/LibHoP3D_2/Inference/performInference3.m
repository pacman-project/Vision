% this is to compute coverage of all vertical layers (3,5, etc)
%cluster3Depths = n2Clusters * n2Clusters * numDisps * 3; - last dimension
%values are: depthMin, depthMax, depthAvr

% X - are combinations 
function [] = performInference3(list_depths, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, n2Clusters, n3Clusters, X, partsOut, coverageOut, displacement, abstractionLevel,  ...
               areLayersRequired, outRoot, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
               is_downsampling, dowsample_rate, elPath, meargeThresh, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, thresh, depthStep, cluster3Depths, dataSetNumber)


    
    hDispl = round(displacement/2);  % half displ
    lenDPW = length(elPath);
    lenCombs = size(X, 1);
    
    isTrim = false;
    isX = false;
    isY = false;
    
    table2 = zeros(n2Clusters, 2);  % required to measure similarities
    for i = 1:n2Clusters
        [clusterX, clusterY] = compute2derivatives(i, nClusters);
        table2(i, 1) = clusterX;
        table2(i, 2) = clusterY;
    end
    
    table3 = uint16(zeros(n2Clusters,n2Clusters,n2Clusters));

    for i = 1:n3Clusters
        cur = partsOut(i,:);  % left centre right
        table3(cur(1),cur(2),cur(3)) = i;   % left centre right
    end
    
    table3Abst = table3;  % assume elements refer to themselve
    

    % merge the selected elements with all elements with distance
    % less than meargeThresh

    % distances has size (lenCombs x n3Clusters)
     X_first = Convert3ToFirstLayer(X, lenCombs, table2);
     partsOutFirst = Convert3ToFirstLayer(partsOut, n3Clusters, table2);  

     distances = Isodata_distances(X_first, partsOutFirst, lenCombs, n3Clusters, false, false);

    for j = 1:lenCombs
        % for each X(j,:) we find the closest element from partsOut
        % table

        if table3(X(j,1), X(j,2), X(j,3)) > 1  % part is selected, no abstraction needed
            continue;
        end

        curDistances = distances(j,:);
        curDistances(curDistances == 0) = 100;  % similarity to itself
        minD = min(curDistances);
        if minD > meargeThresh
            continue;
        end
        idx = find(curDistances == minD);

        if length(idx) > 1
            % Assumption: we have to find the most frequent part
            % according to coverageOut
            coverages = coverageOut(idx);               
            mCoverage = max(coverages);
            idxCov = find(coverages == mCoverage);  % index in the small array
            idx = idx(idxCov(1));
        end

        table3Abst(X(j,1), X(j,2), X(j,3)) = table3(partsOut(idx,1), partsOut(idx,2), partsOut(idx,3));
    end 
    
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end
    
    
    parfor i = 1:lenF 
        
        I    = imread(list_depths{i});
        I = I(:,:,1);
        if dataSetNumber == 1 || dataSetNumber == 3
            mask = [];
        elseif dataSetNumber == 2
            mask = imread(list_mask{i});
        end
        
        I = double(I);
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end

        marks2 = imread(list_El{i});
        
        if dataSetNumber == 1 || dataSetNumber == 3 
            % preliminary processing of the image I
            [I, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                    is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);
        elseif dataSetNumber == 2
        end
        
        [rEl, cEl] = size(mask);  % all three images should be of the same size!
        if r~=rEl || c ~= cEl
            disp('ERROR 20');
        end
        [rEl, cEl] = size(marks2);
        if r~=rEl || c ~= cEl
            disp('ERROR 21');
        end
        if ~is_successfull
            disp('ERROR 22');
        end
        
        marks3 = zeros(r,c);
        % Imarks = zeros(r,c,3);
        
        [rows, cols] = find(marks2 > 0);
        nEl = length(rows);
        
        for j = 1: nEl
            
            % check whether it is close to the boundary
            if rows(j) < displacement + 1  || rows(j) > r - displacement  || cols(j) < displacement + 1 || cols(j) > c - displacement
                continue;
            else
                % otherwise try to match something around this object
                
                central = marks2(rows(j), cols(j));
                depthCentral =  I(rows(j), cols(j));
                
                
                % check what are left and right neighbours
                
                lefts =  [marks2(rows(j), cols(j) - displacement), marks2(rows(j), cols(j) - hDispl)];  % , marks2(rows(j), cols(j) - displ - hDispl)
                rights = [marks2(rows(j), cols(j) + displacement), marks2(rows(j), cols(j) + hDispl)];  % , marks2(rows(j), cols(j) + displ + hDispl)
                depthsLeft =  [I(rows(j), cols(j) - displacement),      I(rows(j), cols(j) - hDispl)];
                depthsRight = [I(rows(j), cols(j) + displacement),      I(rows(j), cols(j) + hDispl)];
                
                indsL = find(lefts > 0);
                indsR = find(rights > 0);

                if (isempty(indsL)) || (isempty(indsR))  % nothing can be matched
                    continue;
                end

                lefts = lefts(indsL);
                rights = rights(indsR);                
                depthsLeft = depthsLeft(indsL);
                depthsRight = depthsRight(indsR);
                
                done = false; % local structure if matched to the vocabulary elements
                ii = 1;
                jj = 1;
                
                while (~done && ii <= length(indsL) && jj <= length(indsR))
                    
                    left = lefts(ii);
                    right = rights(jj);
                    dLeft = depthsLeft(ii);
                    dRight = depthsRight(jj);
                    
                    el = [left, central, right];
                    
                    % check both pairs w.r.t relative depth:
                    
                    relDepthLeft = (dLeft - depthCentral)   / depthStep;
                    relDepthRight = (dRight - depthCentral) / depthStep;
                    
                    is_ok = false;
                    
                    if relDepthLeft >= cluster3Depths(central, left, 1, 1)  && relDepthLeft <= cluster3Depths(central, left, 1, 2)
                        is_ok = true;
                    end
                    if relDepthRight >= cluster3Depths(central, right, 2, 1)  && relDepthLeft <= cluster3Depths(central, right, 2, 2)
                        is_ok = true;
                    end
                    
                    if is_ok  % both pairs are valid, try to match a triple
                    
                    curEl = table3Abst(el(1), el(2), el(3));  % all or nodes are already in this table

                        if curEl ~= 0
                            marks3(rows(j), cols(j)) = curEl;
                            done = true;                     
                        end
                    end
                    
                    jj = jj+1; % increment loop variable
                    if jj > length(indsR)
                        jj = 1;
                        ii = ii+1;
                    end
                    
                end
                               
            end 
        end
        
%         marks3 = marks3+ 100;
%         marks3(marks3 == 100) = 0;
% 
%         imtool(I, [min(min(I)), max(max(I))]);
%         imtool(marks2, [0, n2Clusters]);
%         imtool(marks3, [0, n3Clusters + 100]);

        
        % save the image
        curStr = list_El{i};

        fileName = curStr(lenDPW+1:end);
        outFile = [outRoot, fileName];

        ll = strfind(outFile, '/');
        ll = ll(end); % last position
        folderName = outFile(1:ll);
        b = exist(folderName,'dir');

        if b == 0
            mkdir(folderName);
        end

        marks3 = uint16(marks3);
        imwrite(marks3, outFile, 'png');
        
        if mod(i,200) == 0
            i
        end
        
    end 
        
end