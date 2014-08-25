% this is to compute coverage of all layers after 4 (5, 6, ect.)


%cluster3Depths = n2Clusters * n2Clusters * numDisps * 3; - last dimension
%values are: depthMin, depthMax, depthAvr

% X - are combinations 
function [] = performInferenceNext(list_depths, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                nPrevClusters, nCurClusters, X, partsOut, coverageOut, displ3, displ5, displ7, ...
                areLayersRequired, outRoot, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, elPath, meargeThresh, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster5Depths, fileForVisualization4Layer, ...
                dataSetNumber, layerID)


    is_5_layer = areLayersRequired(5);
    is_6_layer = areLayersRequired(6);
    is_7_layer = areLayersRequired(7);
    is_8_layer = areLayersRequired(8);
    is_inhibition_5_layer = isInhibitionRequired(5);
    is_inhibition_6_layer = isInhibitionRequired(6);
    is_inhibition_7_layer = isInhibitionRequired(7);
    is_inhibition_8_layer = isInhibitionRequired(8);
    
    
    smallOffset = 2;  % half displ
    lenDPW = length(elPath);
    lenCombs = size(X, 1);
    
    isTrim = false;
    isX = false;
    isY = false;
    
    
    table5 = zeros(nPrevClusters,nPrevClusters,nPrevClusters);
    
%   table5 = sparse(zeros(1, nPrevClusters^3));

    tic
    for i = 1:nCurClusters
        cur = partsOut(i,:);  % left centre right    or  top centre bottom
        table5(cur(1),cur(2),cur(3)) = i; 
%         table5((partsOut(i,1) - 1) * (nPrevClusters^2) +  (partsOut(i,2) - 1) * nPrevClusters + partsOut(i,3)) = i;  
%         i
    end
    toc;
    
    table5Abst = table5;  % assume elements refer to themselve
    
    % merge the selected elements with all elements with distance
    % less than meargeThresh

    % distances has size (lenCombs x nPrevClusters)
    if layerID == 5
        [XX] = Convert5ToFirstLayer(X, lenCombs, fileForVisualization4Layer);
        [partsOutXX] = Convert5ToFirstLayer(partsOut, nCurClusters, fileForVisualization4Layer);
    elseif layerID == 6
        [XX] = Convert6ToFirstLayer(X, lenCombs, fileForVisualization4Layer);
        [partsOutXX] = Convert6ToFirstLayer(partsOut, nCurClusters, fileForVisualization4Layer);
    end

    distances = Isodata_distances(XX, partsOutXX, lenCombs, nCurClusters, false, false);

    for j = 1:lenCombs
        % for each X(j,:) we find the closest element from partsOut
        % table

        if table5(X(j,1), X(j,2), X(j,3)) > 1  % part is selected, no abstraction needed
%         if table5((X(j, 1) - 1) * (nPrevClusters^2) +  (X(j, 2) - 1) * nPrevClusters + X(j, 3)) + 0 > 1
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
            % according to coverageOut:   WHICH DOES NOT MAKE TOO MUCH
            % SENCE!
            coverages = coverageOut(idx);               
            mCoverage = max(coverages);
            idxCov = find(coverages == mCoverage);  % index in the small array
            idx = idx(idxCov(1));
        end

        % element X(j,:) should refer to element partsOut(idx)
        table5Abst(X(j,1), X(j,2), X(j,3)) = table5(partsOut(idx,1), partsOut(idx,2), partsOut(idx,3));
        
%         table5Abst((X(j, 1) - 1) * (nPrevClusters^2) +  (X(j, 2) - 1) * nPrevClusters + X(j, 3)) = ...
%             0 + table5((partsOut(idx,1) - 1) * (nPrevClusters^2) +  (partsOut(idx,2) - 1) * nPrevClusters + partsOut(idx,3));
    end 

    
    if dataSetNumber == 1 || dataSetNumber == 3
        list_mask = zeros(1, lenF);
    end
    
    
    for i = 1:lenF   
        
        I = imread(list_depths{i});
        I = I(:,:,1);
        
        if dataSetNumber == 2
            mask = imread(list_mask{i});
        else
            mask = [];
        end
        
        I = double(I);
        
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end

        marks4 = imread(list_El{i});
        
%         marks44 = marks4 + 100;
%         marks44(marks44 == 100) = 0;
%         imtool(marks44, [0, nCurClusters + 100]);
        
        % preliminary processing of the image I
        [I, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                    is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);
                                                                
        [rEl, cEl] = size(mask);  % all three images should be of the same size!
        if r~=rEl || c ~= cEl
            disp('ERROR');
        end
        [rEl, cEl] = size(marks4);
        if r~=rEl || c ~= cEl
            disp('ERROR');
        end
        if ~is_successfull
            disp('ERROR');
        end
        
        marks5 = zeros(r,c);
        
        [rows, cols] = find(marks4 > 0);
        nEl = length(rows);
        
        for j = 1: nEl

            
            % TO BE CHANGED!!!
            
            % check whether it is close to the boundary
            if rows(j) < displ5 + smallOffset + 1  || rows(j) > r - displ5 - smallOffset  || cols(j) < displ5 + smallOffset + 1 || cols(j) > c - displ5 - smallOffset
                continue;
            else
                % otherwise try to match something around this object
                
                central = marks4(rows(j), cols(j));
                depthCentral =  I(rows(j), cols(j));
                
                
                % check what are left and right neighbours
                
                [indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, cols(j), rows(j), displ7, displ5, displ3, smallOffset);
                
                % these are two examples of usage of these function:
                % sub2ind(size(A), [3 2 3 1 2], [3 4 1 3 4], [2 1 2 2 1]);
                % sub2ind(matrixSize, rowSub, colSub);
                
                
                leftInds  = sub2ind(size(marks4), indsYLeft,  indsXLeft);
                rightInds = sub2ind(size(marks4), indsYRight, indsXRight);
                
                lefts       = marks4(leftInds);
                rights      = marks4(rightInds);
                depthsLeft  = I(leftInds);
                depthsRight = I(rightInds);
                
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
                    
                    if relDepthLeft >= cluster5Depths(central, left, 1, 1)  && relDepthLeft <= cluster5Depths(central, left, 1, 2)
                        is_ok = true;
                    end
                    if relDepthRight >= cluster5Depths(central, right, 2, 1)  && relDepthLeft <= cluster5Depths(central, right, 2, 2)
                        is_ok = true;
                    end
                    
                    if is_ok  % both pairs are valid, try to match a triple
                    
                    curEl = table5Abst(el(1), el(2), el(3));  % all or nodes are already in this table
%                   curEl = table5Abst((el(1) - 1) * (nPrevClusters^2) +  (el(2) - 1) * nPrevClusters + el(3));

                        if curEl ~= 0
                            marks5(rows(j), cols(j)) = curEl;
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
        
        marks5 = marks5 + 100;
        marks5(marks5 == 100) = 0;
        imtool(marks5, [0, nCurClusters + 100]);
        
        
%         imtool(I, [min(min(I)), max(max(I))]);
%         imtool(marks5, [0, nPrevClusters + 100]);

        
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

        marks5 = uint16(marks5);
        imwrite(marks5, outFile, 'png');
        
        if mod(i,20) == 0
            i
        end
        
   end
          

        
end