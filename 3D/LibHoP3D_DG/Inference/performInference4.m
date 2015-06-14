% this is to compute coverage of all vertical layers (3,5, etc)
%cluster3Depths = n2Clusters * n2Clusters * numDisps * 3; - last dimension
%values are: depthMin, depthMax, depthAvr

% X - are combinations 
function [] = performInference4(list_depths, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel,...
                n3Clusters, n4Clusters, X, partsOut, coverageOut, displacement,  ...
                outRoot, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, elPath, meargeThresh, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster4Depths, fileForVisualization3Layer, dataSetNumber)

    
    smallOffset = 2;
    lenDPW = length(elPath);
    lenCombs = size(X, 1);
    
    isTrim = false;
    isX = false;
    isY = false;
    
%     table2 = zeros(n2Clusters, 2);  % required to measure similarities
%     for i = 1:n2Clusters
%         [clusterX, clusterY] = compute2derivatives(i, nClusters);
%         table2(i, 1) = clusterX;
%         table2(i, 2) = clusterY;
%     end
    
    table4 = uint16(zeros(n3Clusters,n3Clusters,n3Clusters));

    for i = 1:n4Clusters
        cur = partsOut(i,:);  % left centre right
        table4(cur(1),cur(2),cur(3)) = i;   % left centre right
    end
    
    table4Abst = table4;  % assume elements refer to themselve
    
    
    % merge the selected elements with all elements with distance
    % less than meargeThresh

    % distances has size (lenCombs x n3Clusters)
    [XX] = Convert4ToFirstLayer(X, lenCombs, fileForVisualization3Layer);
    [partsOutXX] = Convert4ToFirstLayer(partsOut, n4Clusters, fileForVisualization3Layer);

    distances = Isodata_distances(XX, partsOutXX, lenCombs, n4Clusters, false, false);

    for j = 1:lenCombs
        % for each X(j,:) we find the closest element from partsOut
        % table

        if table4(X(j,1), X(j,2), X(j,3)) > 1  % part is selected, no abstraction needed
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

        % element X(j,:) should refer to element partsOut(idx)

        table4Abst(X(j,1), X(j,2), X(j,3)) = table4(partsOut(idx,1), partsOut(idx,2), partsOut(idx,3));
    end 
    
    if dataSetNumber == 1 || dataSetNumber == 3  % without it parfor loo de not work :(
        list_mask = zeros(1, lenF);
    end
    
    
    parfor i = 1:lenF 
        
        I    = imread(list_depths{i});
        I = I(:,:,1);
        I = double(I);
        
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end

        marks3 = imread(list_El{i});
        
        if dataSetNumber == 1 || dataSetNumber == 3     % Aim@Shape dataset
            [I, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discRadius, isX, isY, ...
                                isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, [], []);
        elseif dataSetNumber == 2                       % Washington data set
            mask = imread(list_mask{i});
            [I, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY,...
                                isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, [], []);          
        end                                                         
                                                                
        [rEl, cEl] = size(mask);  % all three images should be of the same size!
        if r~=rEl || c ~= cEl
            disp('ERROR 20');
        end
        [rEl, cEl] = size(marks3);
        if r~=rEl || c ~= cEl
            disp('ERROR 21');
        end
        if ~is_successfull
            disp('ERROR 22');
        end
        
        marks4 = zeros(r,c);
        
        [rows, cols] = find(marks3 > 0);
        nEl = length(rows);
        
        for j = 1: nEl
            
            % check whether it is close to the boundary
            if rows(j) < displacement + 1  || rows(j) > r - displacement  || cols(j) < displacement + 1 || cols(j) > c - displacement
                continue;
            else
                % otherwise try to match something around this object
                
                central = marks3(rows(j), cols(j));
                depthCentral =  I(rows(j), cols(j));
                
                
                % check what are bottom and top neighbours (notated as left and right)
                
                lefts =  [marks3(rows(j) - displacement, cols(j)), marks3(rows(j) - displacement, cols(j) - smallOffset), marks3(rows(j) - displacement, cols(j) + smallOffset)];
                rights = [marks3(rows(j) + displacement, cols(j)), marks3(rows(j) + displacement, cols(j) - smallOffset), marks3(rows(j) + displacement, cols(j) + smallOffset)];
                depthsLeft =  [I(rows(j) - displacement, cols(j)),      I(rows(j) - displacement, cols(j) - smallOffset),      I(rows(j) - displacement, cols(j) + smallOffset)];
                depthsRight = [I(rows(j) + displacement, cols(j)),      I(rows(j) + displacement, cols(j) - smallOffset),      I(rows(j) + displacement, cols(j) + smallOffset)];
                
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
                    
                    if relDepthLeft >= cluster4Depths(central, left, 1, 1)  && relDepthLeft <= cluster4Depths(central, left, 1, 2)
                        is_ok = true;
                    end
                    if relDepthRight >= cluster4Depths(central, right, 2, 1)  && relDepthLeft <= cluster4Depths(central, right, 2, 2)
                        is_ok = true;
                    end
                    
                    if is_ok  % both pairs are valid, try to match a triple
                    
                    curEl = table4Abst(el(1), el(2), el(3));  % all or nodes are already in this table

                        if curEl ~= 0
                            marks4(rows(j), cols(j)) = curEl;
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
        
%         marks4 = marks4 + 100;
%         marks4(marks4 == 100) = 0;
%         imtool(marks4, [0, n4Clusters + 100]);
        
        
%         imtool(I, [min(min(I)), max(max(I))]);
%         imtool(marks4, [0, n3Clusters + 100]);

        
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

        marks4 = uint16(marks4);
        imwrite(marks4, outFile, 'png');
        
        if mod(i,20) == 0
            i
        end
        
    end  

        
end