% this is to compute coverage of all layers starting from the first one
%cluster3Depths = n2Clusters * n2Clusters * numDisps * 3; - last dimension
%values are: depthMin, depthMax, depthAvr
% partsOut - selected parts,
% X - combinations observed in the training data!

% X - are combinations 
function [] = performInferenceNext(list_depths, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                nPrevClusters, nCurClusters, nClusters, X, partsOut, partEntropy, displ3, displ5, displ7, ...
                outRoot, isInhibitionRequired, elementRadius, elementType,  ...
                is_downsampling, dowsample_rate, elPrevPath, meargeThresh, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, clusterCurDepths, fileForVisualizationPrevLayer, ...
                dataSetNumber, layerID, cluster1Centres, fieldSize, cluster1Bounds, numSimilar)
    
    emptyCellID = nPrevClusters + 1;  % to tacle empty cells
    
    lenDPW = length(elPrevPath);
    lenCombs = size(X, 1);
    n2Clusters = nClusters ^2;
    
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    fieldCenter = halfFieldSize + 1;
    isTrim = false;
    isX =true;
    isY = true;
    isX_FB = false;
    
    if nPrevClusters > 700
        is_sparse = true;
    else
        is_sparse = false;
    end
    
    if is_sparse
        parfor i = 1:nPrevClusters+1
            tablePrev{i} = sparse(zeros(nPrevClusters+1, nPrevClusters+1));
        end
    else
        tablePrev = zeros(nPrevClusters+1, nPrevClusters+1, nPrevClusters+1);
    end

    
    tic
    for i = 1:nCurClusters
        cur = partsOut(i,:);  % left centre right    or  top centre bottom
        if ~is_sparse
            tablePrev(cur(1), cur(2), cur(3)) = i;
        else
            tablePrev{cur(1)}(cur(2), cur(3)) = i;
        end

    end
    toc;
    tablePrevAbst = tablePrev;  % assume elements refer to themselve
    
    % merge the selected elements with all elements with distance
    % less than meargeThresh

    % distances has size (lenCombs x nPrevClusters)
    downsamplingScheme = 3;
    
    % compute surface descriptors
    [XX, ~, emptyIndicator]            = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer,  ...
                       displ3, displ5, displ7, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths); 
    [partsOutXX, ~, emptyIndicatorPOX] = convertToSurfaceDescriptor(partsOut, nCurClusters, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer,  ...
                       displ3, displ5, displ7, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths);  % this are a first layer descriptors  
                   
    XX(emptyIndicator == 1) = fieldSize(3)*3*depthStep;
    partsOutXX(emptyIndicatorPOX == 1) = fieldSize(3)*3*depthStep;
                  
%   compute distance between surfaces
%   distances = Isodata_distances(XX, partsOutXX, lenCombs, nCurClusters, false, false);
    distances = Integral_distances(XX, partsOutXX, lenCombs, nCurClusters, false, false); 
    counter = 0;
    
    for j = 1:lenCombs
        % for each X(j,:) we find the closest element from partsOut
        % table
        if ~is_sparse
            if tablePrev(X(j,1), X(j,2), X(j,3)) > 1  % part is selected, no need to find closest element from the dictionary       
                continue;
            end
        else
            if tablePrev{X(j,1)}(X(j,2), X(j,3)) > 1
                continue;
            end
        end

        curDistances = distances(j,:);
        curDistances(curDistances == 0) = 100;  % similarity to itself
        minD = min(curDistances);
        idx = find(curDistances == minD);
        
        curCDistances = distances(:,idx);
        % take 5 closets compositions
        [smallestDists, smIDs] = defineSmallestDistances(curCDistances, numSimilar);
        meargeThresh = smallestDists(end) * 4.0;
        
        if minD > meargeThresh
            continue;
        end

        if length(idx) > 1
            % Assumption: we have to find the most frequent part
            % according to coverageOut:   WHICH DOES NOT MAKE TOO MUCH SENSE!
            
%             coverages = coverageOut(idx);               
%             mCoverage = max(coverages);
%             idxCov = find(coverages == mCoverage);  % index in the small array
%             idx = idx(idxCov(1));

            % Assumption 2: take parts with the lowest entropy
            entropies = partEntropy(idx);               
            mCoverage = min(entropies);
            idxCov = find(entropies == mCoverage);  % index in the small array
            idx = idx(idxCov(1));

        end

        % element X(j,:) should refer to element partsOut(idx)
        if ~is_sparse
            tablePrevAbst(X(j,1), X(j,2), X(j,3)) = tablePrev(partsOut(idx,1), partsOut(idx,2), partsOut(idx,3));
            counter = counter + 1;
        else
            tablePrevAbst{X(j,1)}(X(j,2), X(j,3)) = tablePrev{partsOut(idx,1)}(partsOut(idx,2), partsOut(idx,3));
        end
    end 
    
    if dataSetNumber == 1 || dataSetNumber == 3
        list_mask = zeros(1, lenF);
    end
    
    indSS = randperm(lenF);
    
    
    parfor i = 1:lenF 
        
        % save the image
        curStr = list_El{indSS(i)};

        fileName = curStr(lenDPW+1:end);
        outFile = [outRoot, fileName];

        ll = strfind(outFile, '/');
        ll = ll(end); % last position
        folderName = outFile(1:ll);
        
  
        b = exist(folderName,'dir');

        if b == 0
            mkdir(folderName);
        end
        
        b = exist(outFile, 'file');

      if ~b
            
            
        I = imread(list_depths{indSS(i)});
        I = I(:,:,1);
        
        if dataSetNumber == 2
            mask = imread(list_mask{indSS(i)});
        else
            mask = [];
        end
        
        I = double(I);
        
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end

        marksPrev = imread(list_El{indSS(i)});
        
%         marksPrevD = marksPrev + 100;
%         marksPrevD(marksPrevD == 100) = 0;
%         imtool(marksPrevD, [0, nPrevClusters + 100]);
        
        % preliminary processing of the image I
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                    is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, [], [], [], []);
                                                                                                                            
        Ix = Ix(:,:,1);
        Iy = Iy(:,:,1);

        % extend a mask! This is done to tteat areas with empty cells
        mask = extendMaskWithDerivatives(mask, cluster1Bounds, Ix, Iy);
%---------------------------------------------------------                                                               
        
                                                                
        [rEl, cEl] = size(mask);  % all three images should be of the same size!
        if r~=rEl || c ~= cEl
            disp('ERROR 21');
        end
        [rEl, cEl] = size(marksPrev);
        if r~=rEl || c ~= cEl
            disp('ERROR 22');
        end
        if ~is_successfull
            disp('ERROR 23');
        end
        
        marksCur = zeros(r,c);
        
        [rows, cols] = find(marksPrev > 0);
        nEl = length(rows);
        
        % extract a small window arownd each offset
        [indsXOut, indsYOut] = getDispAbs(elementType, elementRadius);
        
        for j=1:nEl
           
                central = marksPrev(rows(j), cols(j));
                depthCentral =  I(rows(j), cols(j));
                
                % check what are left and right neighbours    
                [indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, cols(j), rows(j), displ3, indsXOut, indsYOut); 
                
                % make shure indexes are not out of the image borders
                [indsXLeft, indsYLeft, indsXRight, indsYRight] = checkImageBoundaries(indsXLeft, indsYLeft, indsXRight, indsYRight, r, c);
                
                % NOTE: use a function TestBorder to test two above functions
                
                % these are two examples of usage of these function:
                % sub2ind(size(A), [3 2 3 1 2], [3 4 1 3 4], [2 1 2 2 1]);
                % sub2ind(matrixSize, rowSub, colSub);                
                leftInds  = sub2ind(size(marksPrev), indsYLeft,  indsXLeft);
                rightInds = sub2ind(size(marksPrev), indsYRight, indsXRight);
                
                lefts       = marksPrev(leftInds);
                rights      = marksPrev(rightInds);
                depthsLeft  = I(leftInds);
                depthsRight = I(rightInds);
                leftsMask = mask(leftInds);
                rightsMask = mask((rightInds));
                
                % check is something is empty
                lenEmpLeft = length(leftsMask(leftsMask == 0))/length(indsYLeft);
                lenEmpRight = length(rightsMask(rightsMask == 0))/length(indsYRight);
                
                
                if lenEmpLeft >= 0.5
                    lefts = [lefts, emptyCellID];  % add an empty cell to the list of hypotheses to try
                    depthsLeft = [depthsLeft, depthCentral];
                                       
                end
                if lenEmpRight >= 0.5
                    rights = [rights, emptyCellID];  % add an empty cell to the list of hypotheses to try
                    depthsRight = [depthsRight, depthCentral];
                    
                end
                
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
                    
                    if relDepthLeft >= clusterCurDepths(central, left, 1, 1)  && relDepthLeft <= clusterCurDepths(central, left, 1, 2)
                        is_ok = true;
                    end
                    if relDepthRight >= clusterCurDepths(central, right, 2, 1)  && relDepthLeft <= clusterCurDepths(central, right, 2, 2)
                        is_ok = true;
                    end
                    
                    if is_ok  % both pairs are valid, try to match a triple
                    
                        if ~is_sparse 
                            curEl = tablePrevAbst(el(1), el(2), el(3));  % all OR nodes are already in this table
                        else
                            curEl = tablePrevAbst{el(1)}(el(2), el(3));  % all OR nodes are already in this table
                        end

                        if curEl ~= 0
                            marksCur(rows(j), cols(j)) = curEl;
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
        
 
        
        
%         marksCurD = marksCur + 100;
%         marksCurD(marksCurD == 100) = 0;
%         imtool(marksCurD, [0, nCurClusters + 100]);
%         
%         
%         imtool(I, [min(min(I)), max(max(I))]);
%         imtool(marksCur, [0, nPrevClusters + 100]);

        

        marksCur = uint16(marksCur);
        imwrite(marksCur, outFile, 'png');
        
        if mod(i,2) == 0
            i
        end
        
      end
        
   end
          

        
end



%     is_5_layer = areLayersRequired(5);
%     is_6_layer = areLayersRequired(6);
%     is_7_layer = areLayersRequired(7);
%     is_8_layer = areLayersRequired(8);
%     is_inhibition_5_layer = isInhibitionRequired(5);
%     is_inhibition_6_layer = isInhibitionRequired(6);
%     is_inhibition_7_layer = isInhibitionRequired(7);
%     is_inhibition_8_layer = isInhibitionRequired(8);