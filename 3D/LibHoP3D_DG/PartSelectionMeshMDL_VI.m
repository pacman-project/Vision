% this function performs MDL-based part selection on meshes

% description of the datastructures used here

function [partsOut, coverageOut, lenOut, ORTable] = PartSelectionMeshMDL_VI(list_els, lenFiles, nPrevClusters, iterations, lenSelected, D, layerID)

% list_input, list_els, layerID, nClusters, n2Clusters, ...
%                                     fileForVisualizationLayer, offsetsConventional, depthStep, fieldSize, cluster1Centres, numSimilar, lenSelected)

    % output variables
    coverageOut = zeros(iterations, 1);
    partsOut = zeros(iterations,3);
    lenOut = 0; 
    stat = {};
    
    lenSelected = lenSelected * lenFiles;
                                
    % load files with statistics
    for i = 1:lenFiles 
        strOut = ['Temp/outList_',num2str(i), '.mat'];
        aa = load(strOut);
        stat{i} = aa.outList;
        
        str = ['Temp/outputCoordsAll_A', num2str(i) ,'.mat'];
        aa = load(str);
        outputCoords{i} = aa.outputCoordsAll_A;
%         outputCoords = outputCoords{1}';
    end
    
    stat = [stat{:}]';
    outputCoords = [outputCoords{:}]';
    
    ids1 = stat(:, 1) > 0;
    ids3 = stat(:, 3) > 0;
    idsALL = ids1 & ids3;
    stat = stat(idsALL, :);
    outputCoords = outputCoords(idsALL, :);
    
    strOut = 'Temp/Aggregated.mat';
    aa = load(strOut);
    X = aa.X; 
    lenCombs = size(X, 1);
    
    quantLI{3} = 0.1;
    quantLI{4} = 0.05;
    quantLI{5} = 0.1;
    quantLI{6} = 0.12;
    
    dd = D(:);
    dd = dd(dd > 0);
    quantThresh = quantile(dd, quantLI{layerID});
    
    coverage = zeros(lenCombs, 1);      % just something to initialize
    ORTable = zeros(lenCombs, 1);       % just something to initialize
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
    
    
    %% create a special data Structure that represents all pixels    
    
    numPointsI = zeros(1, lenFiles);  % number of points in each file
    structure = {};
    for i = 1:lenFiles
        fileName1 = list_els{i};
%         if layerID <= 6
            sc = 1;
            strScale = ['scale_', num2str(sc)];
            outFile1 =  [fileName1(1:end-4), strScale, 'NP.mat'];
%         else
%             outFile1 =  [fileName1(1:end-4), 'NP.mat'];
%         end
        str = ['layer', num2str(layerID-1)];
        outFile1 = strrep(outFile1, str ,'layer2');
        aa = load(outFile1);  %'numPoints'
        numPoints = aa.numPoints;
        curLength = sum(numPoints);
        numPointsI(i) = curLength;
        aa = zeros(curLength, 1);
        structure{i} = aa;
    end
    
    
    %% read as much as possible of the data from files to the structures
    for i = 1:lenFiles
        fileStructure{i} = {};
        fileNamesStructure{i} = {};
    end
    for i = 1:lenFiles
        fileName2 = list_els{i};
        if layerID == 3
            sc = 1;
            strScale = ['scale_', num2str(sc)];
            outFile2 = [fileName2(1:end-4), strScale, 'PS.mat'];
        else
            outFile2 = [fileName2(1:end-4), 'PS.mat'];
        end
        fileNamesStructure{i} = outFile2;
    end
    if lenFiles < 12
        for i = 1:lenFiles
            dd = load(fileNamesStructure{i});
            fileStructure{i} = dd.NeiOut;
        end
    end
    
    %%
    [coverage] = recomputeCoverageMeshes(stat, outputCoords, coverage, X, 1, ...
                                        list_els, lenFiles, nPrevClusters + 1, lenCombs, is_sparse, triples, lenCombs, numPointsI, layerID, ORTable);
    save('Temp/covFirst.mat', 'coverage');                                
    dd = load('Temp/covFirst.mat');
    coverage = dd.coverage;
                                    
    % sort according to coverage
    [coverage, indss] = sort(coverage, 'descend');
    X = X(indss, :);   
    X = int16(X);
    D = D(indss, indss);

    iterations = min(iterations, size(X, 1));

    for k = 1:iterations   %kk:iterations+kk   % for each part 
        
        str = ['element - ', num2str(k)];
        disp(str);
        str = ['coverage - ', num2str(coverage(k))];
        disp(str);
        
        if coverage(k) == 0
            return;  % all parts are selected
        end 
        
        
        %% merge the selected elements with all elements with distance less than meargeThresh
        
%         distance = Integral_distances(X_first, X_first(k,:), lenCombs, 1, false, false);
%         [smallestDists, smIDs] = defineSmallestDistances(distance, numSimilar);
%         curPart = X(smIDs, :); 
%         numSimilarUpd = size(curPart, 1); % may be a bit different from the previous one
         
        curPart = X(k, :);
        numSimilarUpd = 1;
        
        % find position of curParts in all images
%         try
%             inds = arrayfun(@isMemberMy, stat(:,1), stat(:,2), stat(:,3));
%             elPositions = outputCoords(inds, :);
%             
%         catch exception
            
        inds = ismember(stat, curPart, 'rows');
        elPositions = outputCoords(inds, :);
        
%         end


        % elPositions  - position of this element in the models: [modelID, centreIDx, leftIDx, rightIDx]

        % project this part to the images
        % curElPosition is already sorted by image number
        firstColumn = elPositions(:,1);
        disp('Projecting the selected part to training data');
        

        for jj = 1:lenFiles
            
            % extract all related to the image
            % 1) find the current position in the elPositions
            ind = find(firstColumn == jj);
            lenE = length(ind);
            
            if lenE > 0
                curStruct = structure{jj};
                
                % read the file
                
                if isempty(fileStructure{jj})
                    outFile2 = fileNamesStructure{jj};
                    dd = load(outFile2);
                    NeiOut1 = dd.NeiOut;
                else
                    NeiOut1 = fileStructure{jj};
                end

                % the right model is now open
                temp = elPositions(ind,2:4);
                aa = NeiOut1(temp(:));
                idsss = [aa{:}];
                dd = length(idsss); 
                step = 10^7;
                numSteps = length(1:step:dd);
                curStep = 0;
                for it = 1:step:dd
                    curStep = curStep+1;
                    beginIT = it;
                    endIT = min(it+step, dd);
                    curStruct(idsss(beginIT:endIT)) = 1;
                end
          
                structure{jj} = curStruct;
            end
        end
        
        coverageOut(k) = coverage(k);
        partsOut(k, :) = X(k,:);
        lenOut = lenOut + 1;
        % check if there are some parts that are can be combined with that
        % one by OR nodes
        idsOR = find(D(k, :) <= quantThresh);
        idsOR = idsOR(ORTable(idsOR) == 0);
        ORTable(idsOR) = k;
 
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
                                        list_els, lenFiles, nPrevClusters + 1, lenRec, is_sparse, triples, lenCombs, numPointsI, layerID, ORTable);
                                        
        % sort parts according to their coverage (sort only the rest of the array)
        [coverage(k+1:end), indsC] = sort(coverage(k+1:end), 'descend');
        indsC = indsC+k;
        X(k+1:end,:) = X(indsC, :);
        ORTable(k+1:end) = ORTable(indsC);
        dTemp = D(k+1:end,k+1:end);
        D(k+1:end,k+1:end) = dTemp(indsC-k, indsC-k);
        
        
%         X_first(k+1:end,:) = X_first(indsC, :);
        
        if k ~= size(X, 1)
            if coverage(k+1)<= lenSelected  % maximal remaining coverage is less than we expect
                idsOT = ORTable>0;
                partsOut = X(idsOT, :);
                coverageOut = coverage(idsOT);
                lenOut = nnz(idsOT);
                ORTable = ORTable(idsOT);
                return;
            end
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
                                        list_els, lenF, nPrevClusters, numRecompute, is_sparse, triples, lenCombs, numPointsI, layerID, ORTable)
    
    disp('recomputing parts coverage...');
    
    X = uint32(X);
    statistics = uint32(statistics);
    
    lenStat = size(statistics, 1);
    
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
    for ii = 1:lenF    % PARFOR
        
        %% read the model
        if isempty(fileStructure{ii})
            outFile3 = fileNamesStructure{ii};
            ddd = load(outFile3);
            NeiOut = ddd.NeiOut;
        else
            NeiOut = fileStructure{ii};
        end
        
        curStr = structure{ii};
        localSumCov = zeros(numRecompute, 1);
        
        % extract parts corresponding to this image from outputCoords
        indsI = secondColumn == ii;
        els = outputCoords(indsI, :);  % all elements of this Model
        firstCol = els(:, 1);
        [~, inddds] = sort(firstCol); % sort according to the elementID
        els = els(inddds, :);
        lenEls = size(els, 1);  % [elID, modelID, CentralIDx, LeftIDx, RightIDx] 
        
        
        if lenEls > 0   % for some images there are no detected elements   
            prevEl = els(1,1);
            jPrev = 1;
            tempStr = zeros(size(curStr));
            
            for j = 1:lenEls
                
                curEl = els(j,1);

                if (curEl ~= prevEl) || j == lenEls % evaluate coverage of the prevEl
                    elsTemp = els(jPrev:j-1,3:5);
                    ids = [NeiOut{elsTemp(:)}]';
                    tempStr(ids) = 1;
                    cover = sum(tempStr & ~curStr);
                    localSumCov(prevEl - startX + 1) = localSumCov(prevEl - startX + 1) + cover;
                    tempStr(:) = 0;
                    prevEl = curEl;
                    jPrev = j;
                end

                if mod(j, 10^6) == 0
                    strDisp = ['File ', num2str(ii),': ',num2str(j), ' out of ', num2str(lenEls)];
                    disp(strDisp);
                end
            end
           
                
        end

        
        sumCov(:, ii) = localSumCov;
    end
    % now write information from coverageTemp to coverage
    coverage(startX:startX+numRecompute-1) = sum(sumCov,2);
    coverage(ORTable > 0) = 0;
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



%         if lenEls > 0   % for some images there are no detected elements   
%             prevEl = els(1,1);
%             
%             for j = 1:lenEls
%                 
%                 curEl = els(j,1);
% 
%                 if curEl ~= prevEl % evaluate coverage of the prevEl
%                     cover = sum(curStr(:,2) & ~curStr(:,3));
%                     localSumCov(prevEl - startX + 1) = localSumCov(prevEl - startX + 1) + cover;
%                     curStr(:,2) = 0;
%                     prevEl = curEl;
%                 end
%                 
%                 ids = [NeiOut{els(j,3)}, NeiOut{els(j,4)}, NeiOut{els(j,5)}]';
%                 curStr(ids, 2) = 1;
%                 
%                 if mod(j, 10^6) == 0
%                     strDisp = [num2str(j), ' out of ', num2str(lenEls)];
%                     disp(strDisp);
%                 end
%             end
%             
%             if j == lenEls
%                 cover = sum(curStr(:,2) & ~curStr(:,3));
%                 localSumCov(curEl - startX + 1) = localSumCov(curEl - startX + 1) + cover;
%             end
%                 
%         end




%         fileName3 = list_els{ii};
%         if layerID == 3 
%             scC = 1;
%             strScaleC = ['scale_', num2str(scC)];
%             outFile3 =  [fileName3(1:end-4), strScaleC, 'PS.mat'];
%         else
%             outFile3 =  [fileName3(1:end-4), 'PS.mat'];
%         end
%         st3 = load(outFile3);
%         NeiOut = st3.NeiOut;













