% this function performs MDL-based part selection on meshes

% description of the datastructures used here

function [partsOut, coverageOut, lenOut] = PartSelectionMeshMDL(statisticsLayer, statisticsLayerSieved, statisticsLayerAggregated, list_input, list_els, layerID, nClusters, n2Clusters, ...
                                    fileForVisualizationLayer, offsetsConventional, depthStep, fieldSize, cluster1Centres, iterations, numSimilar, lenSelected)

    % output variables
    coverageOut = zeros(iterations, 1);
    partsOut = zeros(iterations,3);
    lenOut = 0;                            
                                
%     load(statisticsLayer);
    st1 = load(statisticsLayerSieved);      % 'statistics', 'clusterCurDepths', 'outputCoords', 'outputFrames', '-v7.3');
    st2 = load(statisticsLayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
    
    statistics = st1.statistics;
    outputCoords = st1.outputCoords;
    clusterCurDepths = st1.clusterCurDepths;
    X = st2.X;
    clear('st1', 'st2');
    
    halfFieldSize = fieldSize/2;
    fieldCenter = halfFieldSize;
    
    lenCombs = size(X, 1);
    lenFiles = length(list_input);
   
    downsamplingScheme = 3;
    
    [X_first, nPrevClusters, emptyIndicator] = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationLayer{layerID-1},  ...
                       offsetsConventional, depthStep, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths);  % this is a first layer descriptor
    
    
    X_first(emptyIndicator == 1) = fieldSize(3)*3*depthStep;
    coverage = zeros(lenCombs, 1);      % just something to initialize
    coverage = double(coverage);
    
    % this is done to speed up a function of reconstruction
    
    if nPrevClusters > 800
        is_sparse = true;
    else
        is_sparse = false;
    end
    
    % create a matrix once
    if is_sparse
        parfor i = 1:nPrevClusters
            triples{i} = sparse(zeros(nPrevClusters + 1, nPrevClusters + 1));
        end
    else
        triples = [];  % otherwise create a matrix in the function recomputeCoverage
    end 
    
    indsLine = [2,1,4];
    stat = statistics(:, indsLine);   %  stat becomes: [left, central, right]
    
    %% create a special data Structure that represents all pixels    
    
    numPointsI = zeros(1, lenFiles);  % number of points in each file
    structure = {};
    parfor i = 1:lenFiles
        fileName = list_els{i};
        outFileAP = [fileName(1:end-4), 'AP.mat']; 
        aa = load(outFileAP);  % 'VAll', 'NAll', 'numPoints', 'areCentral' % THESE ARE ADDITIONAL POINTs generated from the mesh
        numPoints = aa.numPoints;
        curLength = sum(numPoints);
        numPointsI(i) = curLength;
        aa = zeros(curLength, 3);
        aa(:,1) = i;
        structure{i} = aa;
    end
    
    [coverage] = recomputeCoverageMeshes(stat, outputCoords, coverage, X, 1, ...
                                        list_input, list_els, lenFiles, nPrevClusters + 1, lenCombs, is_sparse, triples, lenCombs, numPointsI);
                                    
    % sort according to coverage
    [coverage, indss] = sort(coverage, 'descend');
    X = X(indss, :);
    X_first = X_first(indss, :);   
    X = int16(X);
%     stat = gpuArray(stat);
    
    
    for k = 1:iterations   %kk:iterations+kk   % for each part 
        
        str = ['element - ', num2str(k)];
        disp(str);
        str = ['coverage - ', num2str(coverage(k))];
        disp(str);
        
        if coverage(k) == 0
            return;  % all parts are selected
        end 
        
        
        % merge the selected elements with all elements with distance
        % less than meargeThresh
        
        distance = Integral_distances(X_first, X_first(k,:), lenCombs, 1, false, false);
        [smallestDists, smIDs] = defineSmallestDistances(distance, numSimilar);
        
        
        curPart = X(smIDs, :); 
        numSimilarUpd = size(curPart, 1); % may be a bit different from the previous one
        
        % find position of curParts in all images
        try
            inds = arrayfun(@isMemberMy, stat(:,1), stat(:,2), stat(:,3));
            elPositions = outputCoords(inds, :);
            
        catch exception
            
            inds = ismember(stat, curPart, 'rows');
            elPositions = outputCoords(inds, :);
        end
        
        % elPositions  - position of this element in the models: [modelID, centreIDx, leftIDx, rightIDx]

        % project this part to the images
        % curElPosition is already sorted by image number
        firstColumn = elPositions(:,1);
        disp('Projecting the selected part to training data');
        
        parfor jj = 1:lenFiles
            
            % extract all related to the image
            % 1) find the current position in the elPositions
            ind = find(firstColumn == jj);
            lenE = length(ind);
            
            if lenE > 0
                curStruct = structure{jj};
                
                % read the file
                fileName = list_els{jj};
                outFileMod = [fileName(1:end-4), '.mat'];
                dd = load(outFileMod);   % 'Vout', 'Nout', 'likelihoods', 'darFrames', 'partIDs', 'NeiOut', 'pointIndexing'
                pointIndexing1 = dd.pointIndexing';
                NeiOut1 = dd.NeiOut;

                % the right model is now open
                for kk = 1:lenE  % fill all positions for this image 
                    temp = elPositions(ind(kk),2:4);
                    for kkk = 2:4  % central,left and right elements
                        centralIDxx = temp(1);
                        leftIDxx = temp(2);
                        rightIDxx = temp(3);
                        cc_cent = pointIndexing1(:,centralIDxx);  % fid, pointID
                        cc_left = pointIndexing1(:,leftIDxx);
                        cc_right = pointIndexing1(:,rightIDxx);
                        idsss = [NeiOut1{cc_cent(1)}{cc_cent(2)}, NeiOut1{cc_left(1)}{cc_left(2)}, NeiOut1{cc_right(1)}{cc_right(2)}]';
                        curStruct(idsss, 3) = 1;
                    end
                    
                end
                structure{jj} = curStruct;
            end
        end
        
        coverageOut(k) = coverage(k);
        partsOut(k, :) = X(k,:);
        lenOut = lenOut + 1;
        
 
        %------------------------------------------------------------------
        % update coverage of the following parts
        lenRec = lenCombs - k;
%         fullRecFreq = 3;
%         
%         if mod(k, fullRecFreq) == 0 || (coverage(k+1)<= lenSelected) % full recompute
%             lenRec = lenCombs - k;
%         elseif mod(k, 2) == 0
%             lenRec = 30;
%         end

        [coverage] = recomputeCoverageMeshes(stat, outputCoords, coverage, X, k+1, ...
                                        list_input, list_els, lenFiles, nPrevClusters + 1, lenRec, is_sparse, triples, lenCombs, numPointsI);
                                        
        
        % sort parts according to their coverage (sort only the rest of the array)
        [coverage(k+1:end), indsC] = sort(coverage(k+1:end), 'descend');
        indsC = indsC+k;
        X(k+1:end,:) = X(indsC, :);
        X_first(k+1:end,:) = X_first(indsC, :);
        
        if coverage(k+1)<= lenSelected  % maximal remaining coverage is less than we expect
            return;
        end
            
        save('Temp/partsOut.mat', 'partsOut');
    end
    
%     pointsOverall = 0;
%     pointsCovered = 1;
%     
%     for aa = 1:lenFiles  
%         
%         
%     end
    a = 2;
%     fileName = list_els{i};
%     outFileM = [fileName(1:end-4), '.mat'];
%     outFileAP = [fileName(1:end-4), 'AP.mat']; 
% 
%     [~, F, ~] = meshRead(list_input{i});
%     load(outFileM);   % 'Vout', 'Nout', 'likelihoods', 'darFrames', 'partIDs', 'NeiOut', 'pointIndexing'
%     load(outFileAP);  % 'VAll', 'NAll', 'numPoints', 'areCentral' 
     



%% This is a function to recompute parts coverage

function [coverage] = recomputeCoverageMeshes(statistics, outputCoords, coverage, X, startX, ...
                                        list_input, list_els, lenF, nPrevClusters, numRecompute, is_sparse, triples, lenCombs, numPointsI)
    
    disp('recomputing parts coverage...');
    
    X = uint32(X);
    statistics = uint32(statistics);
    
    lenStat = size(statistics, 1);
    lenCombs = length(coverage);
    
    if startX+numRecompute-1 > lenCombs
        numRecompute = lenCombs - startX + 1;
    end
          
    coverage(startX:startX+numRecompute-1) = coverage(startX:startX+numRecompute-1) * 0;

    % create a matrix
    if isempty(triples)
        triples = zeros(nPrevClusters, nPrevClusters, nPrevClusters);
    end
 
    % create a table of format [el, modelID, centralIdx, LeftIDx, RightIDx]
    outputCoords = [zeros(lenStat, 1), outputCoords];
    
    if is_sparse 
        % fill a table triples
        for t = startX:startX+numRecompute-1  % lenCombs
            triples{X(t, 1)}(X(t,2), X(t,3)) = t;
        end
        
        % fill a table outputCoords
        for t = 1:lenStat
            outputCoords(t,1) = triples{statistics(t,1)}(statistics(t,2), statistics(t,3)); 
        end
               
    else    % if it is NOT sparse
        
        % fill a table triples
        ids = startX:1:startX+numRecompute-1;
        indSS = sub2ind(size(triples), X(ids,1), X(ids,2), X(ids,3)); 
        triples(indSS) = ids;
    
%         for i = startX:startX+numRecompute-1  % lenCombs
%             triples(X(i, 1), X(i,2), X(i,3)) = i;
%         end
        
        % fill a table outputCoords
        
        ids = 1:1:lenStat;
        indSS = sub2ind(size(triples), statistics(ids,1), statistics(ids,2), statistics(ids,3)); 
        outputCoords(ids, 1) = triples(indSS);
        
%         for i = 1:lenStat 
%             outputCoords(i,1) = triples(statistics(i,1), statistics(i,2), statistics(i,3));  
%         end 
    end

    % to exclude already selected parts (and those which do not need to be recomputed)
    firstCol = outputCoords(:,1);
    indEx = firstCol ~= 0;
    outputCoords = outputCoords(indEx, :);

    
    %----------------------------------------------------------------------
    secondColumn = outputCoords(:,2); % [el, Im, x,y]
    [secondColumn, idx] = sort(secondColumn, 'ascend'); % sort to speed up
    outputCoords = outputCoords(idx, :);
    
    disp('computing coverage through all images...');
    sumCov = zeros(numRecompute, lenF);
    
    parfor ii = 1:lenF  % for each model
        
        %% read the model
        fileNameC = list_els{ii};
        outFileM = [fileNameC(1:end-4), '.mat'];
%       outFileAP = [fileNameC(1:end-4), 'AP.mat'];
%       load(outFileAP);  % 'VAll', 'NAll', 'numPoints', 'areCentral'

        curStr = structure{ii};
        localSumCov = zeros(numRecompute, 1);

        st3 = load(outFileM);   % 'Vout', 'Nout', 'likelihoods', 'darFrames', 'partIDs', 'NeiOut', 'pointIndexing'
        partIDs = st3.partIDs;
        NeiOut = st3.NeiOut;
        pointIndexing = st3.pointIndexing';

        
        % extract parts corresponding to this image from outputCoords
        indsI = secondColumn == ii;
        els = outputCoords(indsI, :);  % all elements of this Model
        firstCol = els(:, 1);
        [~, inddds] = sort(firstCol); % sort according to the elementID
        els = els(inddds, :);
        lenEls = size(els, 1);  % [elID, modelID, CentralIDx, LeftIDx, RightIDx] 
        
        if lenEls > 0   % for some images there are no detected elements   
            prevEl = els(1,1);
            
            for j = 1:lenEls
                
                curEl = els(j,1);

                if curEl ~= prevEl % evaluate coverage of the prevEl
                    cover = sum(curStr(:,2) & ~curStr(:,3));
                    localSumCov(prevEl - startX + 1) = localSumCov(prevEl - startX + 1) + cover;
                    curStr(:,2) = 0;
                    prevEl = curEl;
                end
                
                centralIDx = els(j,3);
                leftIDx = els(j,4);
                rightIDx = els(j,5);
                c_cent = pointIndexing(:,centralIDx);  % fid, pointID
                c_left = pointIndexing(:,leftIDx);
                c_right = pointIndexing(:,rightIDx);
                ids = [NeiOut{c_cent(1)}{c_cent(2)}, NeiOut{c_left(1)}{c_left(2)}, NeiOut{c_right(1)}{c_right(2)}]';
                curStr(ids, 2) = 1;
            end
            
            if j == lenEls
                cover = sum(curStr(:,2) & ~curStr(:,3));
                localSumCov(curEl - startX + 1) = localSumCov(curEl - startX + 1) + cover;
            end
                
        end
        
        sumCov(:, ii) = localSumCov;
    end
    % now write information from coverageTemp to coverage
    coverage(startX:startX+numRecompute-1) = sum(sumCov,2);

end

% nested function impemented on gpu
function isTrue = isMemberMy(col1, col2, col3)
    isTrue = false;
    for jjj = 1:numSimilarUpd
        isTrue = isTrue | (curPart(jjj,1) == col1 & curPart(jjj,2) == col2 & curPart(jjj,3) == col3);
    end
end
        
        
    reset(gpuDevice(1));

end














