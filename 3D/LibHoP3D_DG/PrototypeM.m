% % this is prototype of the whole corner and edge learning-inference
    
function [] = PrototypeM()

    dataSetNumber = 3;
    nClusters{1} = 9;
    nClusters{2} = nClusters{1}^2;
    nClusters{3} = 25;
    subsetPercent = 1.0;
    is_subset = false;
    structureTensorRad = 16;
    
    is_inference{1} = 0;
    
    [cluster1Centres, cluster1Bounds, thresh] = defineFirstLayerQuantiles(nClusters{1}, dataSetNumber); 
    [filtOptions, zeroThresh] = loadFilteringParameters(dataSetNumber);
    options = GetOptions(filtOptions);
    options.elementRadius = 3;

    input_path = 'D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll';
    output_path{2} = 'D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\layer2\';
    output_path{3} = 'D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\layer3\';
    output_path{4} = 'D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\layer4\';
    
    outputNormalFolder = 'D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\layer2_Normals\';
    visualizationFolder{3} = 'D:\Visualized vocabulary\Vladislav_STD\3_layer_loc_new\';
    visualizationFolder{4} = 'D:\Visualized vocabulary\Vladislav_STD\3_layer_loc_new\';
    
    if ~exist(output_path{2},'dir')
        mkdir(output_path{2});
    end
    if ~exist(output_path{3},'dir')
        mkdir(output_path{3});
    end
    if ~exist(output_path{4},'dir')
        mkdir(output_path{4});
    end
    if ~exist(outputNormalFolder,'dir')
        mkdir(outputNormalFolder);
    end
    
    zPath = 'D:\LibHoP3D_DG\settings\calibration1.mat';
    load(zPath);

    [list_input, list_mask, ~, lenF] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);
    
    if is_inference{1}
        FirstLyerInference(list_input, lenF, zScale, options, cluster1Bounds, nClusters{1}, output_path, outputNormalFolder);
    end
 
    %% on the next step we perform learning - inference of the 3rd layer

    for layerID = 3:3
        
        [list_els2, ~, ~,lenF] = extractFileListGeneral(output_path{2}(1:end-1), is_subset, subsetPercent, dataSetNumber);
        [list_els, ~, ~, lenF] = extractFileListGeneral(output_path{layerID-1}(1:end-1), is_subset, subsetPercent, dataSetNumber);
        
%         learningInferenceNext(list_input, list_els, list_els2, lenF, zScale, options, output_path{layerID}, outputNormalFolder, ...
%             nClusters{1}, nClusters{layerID - 1}, layerID, structureTensorRad);

%         [list_els, ~, ~, lenF] = extractFileListGeneral(output_path{layerID}(1:end-1), is_subset, subsetPercent, dataSetNumber);
%         [partIDs, nParts] = symmetricTrial(nClusters{2}); 
%         recomputePartIds(list_els, lenF, output_path{3}, nClusters{2}, partIDs, nClusters{1});  % This is to make all part ids of the layer in range from 0 to N
% 
%         [list_elsPrev, ~, ~, ~] = extractFileListGeneral(output_path{2}(1:end-1), is_subset, subsetPercent, dataSetNumber);
%         [list_elsCur, ~, ~, ~] = extractFileListGeneral(output_path{3}(1:end-1), is_subset, subsetPercent, dataSetNumber);
% 
%         % collect co-occurrence statistics for part visualization
%         statAll = getCooccurrenceStat(list_elsPrev, list_elsCur, list_input, zScale, structureTensorRad);
        if layerID == 3
            load('statMirL_3.mat'); % 'statAll'
            load('clusterTable_3.mat');
            lenST = length(statAll);
            relPos = zeros(lenST,3);

            %% visualization of the parts
            for i = 1:lenST
                temp = statAll{i};
                relPos(i, :) = sum(temp, 2)'/size(temp,2);
            end
            visualizeParts(relPos, clusterCentresT, visualizationFolder{3}, cluster1Centres, nClusters{1}, layerID);
        end
    end
    
end

function [] = FirstLyerInference(list_input, lenF, zScale, options, cluster1Bounds, nClusters, output_path, outputNormalFolder)

    gf = fspecial('gaussian',5,1);
    dx = 0.5 * [-1, 0, 1];
    dy = dx';
    elementRadius = options.elementRadius;
    
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    
    parfor i = 1:lenF
        
        I = imread(list_input{i});
%         I = imread('D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\99_1_1_3_2_5.png');
        I = double(I(:,:,1));
        mask = I > 0;
        I = I * zScale;
        I = imfilter(I, gf);
%         Ix = imfilter(I, dx);
%         Iy = imfilter(I, dy);
        [r,c,ch] = size(I);

        outImage = zeros(r,c);
        Normals = zeros(r,c,3);
        

        % zScale
        
        for y = 1:2:r
            for x = 1:2:c
                if mask(y,x) == 0 || I(y,x) < options.minValue
                    continue;
                end
                
                [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, elementRadius, options.is2D, options.minValue, options.maxValue);
                if length(indsXOut) < options.ignoreThresh{1};
                    continue;
                end

                % now we convert these coordinates to the local frame of
                % reference (centered at the point [j,i,curDepth])
                xs = indsXOut - x;
                ys = indsYOut - y;
                zs = depths - I(y,x);
                if length(xs) < options.ignoreThresh{1} + 2
                    continue;
                end
                
                [absErr, Norm, is_ok] = PlanarTestLeastSquares(xs, ys, zs);
                if is_ok
                    angleX = 90 - 180 * acos(dot(Xtemp,Norm)/norm(Norm))   /pi;
                    angleY =  90 - 180 * acos(dot(Ytemp,Norm)/norm(Norm))   /pi;
                    
                    clusterX = define1Cluster(angleX, cluster1Bounds, nClusters);
                    clusterY = define1Cluster(angleY, cluster1Bounds, nClusters);
                    if clusterX == 0 || clusterY == 0
                        continue;
                    end
                    partID = compute2elementIndex(clusterX, clusterY, nClusters);
                    outImage(y,x) = partID;
                    Normals(y,x,1) = Norm(1); Normals(y,x,2) = Norm(2); Normals(y,x,3) = Norm(3);
                end

            end
        end
        % save the part ID and Normals to the files
        [~, fileName] = getFolderName(list_input{i});
        outFileNormal = [outputNormalFolder, fileName];
        outFileNormal(end-3:end) = '.mat';
        outFile = [output_path{2}, fileName];
        parSave(outFileNormal, Normals);
        
        outImage = uint8(outImage);
        imwrite(outImage, outFile, 'png');
        if mod(i,10) == 0
            disp(i);
        end
    end
end

function [tablelParts, tableCoords] = learningInferenceNext(list_input, list_els, list_els2, lenF, zScale, options, output_path, outputNormalFolder, nClusters, nPrevClusters, layerID, structureTensorRad)
    
    tablelParts = zeros(10, 3);
    tableCoords = zeros(10, 9);
    curr = 0;

    step = 2;
    nonVanishThresh = 0.27;
    nonVanishThresh3 = 0.23;
    gf = fspecial('gaussian',5,1);
    dx = 0.5 * [-1, 0, 1];
    dy = dx';

    ignoreThresh = 4; % to prevent creation of pairs on the same side of the egde

    [partIDs, nParts] = symmetricTrial(nPrevClusters);
    offset = 12;
    
    % look up table
    
    if layerID == 3
        LUTable = zeros(2,nPrevClusters);
        for i = 1:nPrevClusters
            [clusterX, clusterY] = compute2derivatives(i, nClusters);
            LUTable(1,i) = clusterX;
            LUTable(2,i) = clusterY;
        end
    end
    
    for i = 9:lenF
        
%         I = imread('D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\99_1_1_3_2_5.png');
%         marksPrev = imread('D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\layer2\99_1_1_3_2_5.png');      
%         [~, fileName] = getFolderName('D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\99_1_1_3_2_5.png');

%         str = ['start ', num2str(i)];
%         disp(str);
        
        I = imread(list_input{i});
        marksPrev = imread(list_els{i});
        if layerID ~= 3
            marks2 = imread(list_els2{i});
        else 
            marks2 = marksPrev;
        end
        
        if max(max(marksPrev)) == 0
            continue;
        end
        
        [~, fileName] = getFolderName(list_input{i});
        NormalFile = [outputNormalFolder, fileName];
        NormalFile(end-3:end) = '.mat';
        Normals = load(NormalFile);
        Normals = Normals.Normals;
        
        outFile = [output_path, fileName];

        I = double(I(:,:,1));
%         marksPrev = imread(list_els{i});
        mask = I > 0;
        I = I * zScale;
        I = imfilter(I, gf);
        Ix = imfilter(I, dx);
        Iy = imfilter(I, dy); 
        [r,c,ch] = size(I);

        outMarks = zeros(r,c);
        %% specify the structure tensor for each location 

        for x =  structureTensorRad+1 : step: c - structureTensorRad - 1
            for y = structureTensorRad+1 : step: r - structureTensorRad - 1
                
                if mask(y,x) == 0
                    continue;
                end 
                
                [~, indsX, indsY, ~, ~] = GetPartsNeighbour(I, marks2, x, y, I(y,x), structureTensorRad, false);
                numPoint = length(indsX);
                if numPoint < 3
                    continue;
                end
                
                NormalsTemp = zeros(numPoint,3);
                for j = 1:numPoint
                    NormalsTemp(j, 1) = Normals(indsY(j), indsX(j), 1);
                    NormalsTemp(j, 2) = Normals(indsY(j), indsX(j), 2);
                    NormalsTemp(j, 3) = Normals(indsY(j), indsX(j), 3);
                end
                
                % compute a structure tensor here
                T = (NormalsTemp' * NormalsTemp);
                T = T/numPoint;
                [V,D] = eig(T);
                b = diag(D);
                bs = abs(b);
                [bs, ids] = sort(bs, 'descend');
                b = b(ids);
                
%                 imD(y,x,1) = abs(b(1))*2; imD(y,x,2) = abs(b(2))*2; imD(y,x,3) = abs(b(3))*2;

                cond1 = (abs(b(1)) >= nonVanishThresh) && abs(b(2)) >= nonVanishThresh && abs(b(3)) <= nonVanishThresh && layerID == 3;
                cond2 = (abs(b(1)) >= nonVanishThresh) && abs(b(2)) >= nonVanishThresh && abs(b(3)) <= nonVanishThresh && layerID == 4;
                cond3 = (abs(b(1)) >= nonVanishThresh3) &&abs(b(2)) >= nonVanishThresh3 && abs(b(3))>= nonVanishThresh3 &&layerID == 4;

                %% edge parts
                if cond1  % this is an edge
%                     % 1 compute principal directions
                    indIm = sub2ind(size(marks2), indsY, indsX);
                    depths = I(indIm);
                    points = [indsX', indsY', depths'];
                    pointCentral = [x,y, I(y,x)];
                    pointsLocal = points - repmat(pointCentral, [size(points, 1),1]);
                    
%                     figure;
%                     scatter3(pointsLocal(:,1), pointsLocal(:,2), pointsLocal(:,3));
%                     hold on
  
                     % 2 use this points to estiamte the Darboux frame
                    tangentX = [1, 0, Ix(y,x)];
                    tangentY = [0, 1, Iy(y,x)];
                    Norm = cross(tangentX, tangentY);
                    Norm = Norm/norm(Norm);
                    [V, values] = computeDarbouxFrame(Norm, pointsLocal(:, 1), pointsLocal(:, 2), pointsLocal(:, 3));
                                   
%                     vectColors = eye(3);
%                     vecLen = 5.0;
%                     plotFrame(V, [0,0,0], vecLen, vectColors, 3);
%                     axis equal;
%                     
                    dir = V(1:2, 2);
                    pointLeft = round([x,y] + dir' * offset);
                    pointRight = round([x,y] - dir' * offset);
                    
%                     if rand < 0.1
%                         I(y,x) = 0;
%                         I(pointLeft(2), pointLeft(1)) = 0;
%                         I(pointRight(2), pointRight(1)) = 0;
%                     end

                    [marksSelL, ~, ~, ~, ~] = GetPartsNeighbour(I, marksPrev, pointLeft(1),  pointLeft(2),  I(pointLeft(2),  pointLeft(1)), ceil(structureTensorRad/3), false);
                    [marksSelR, ~, ~, ~, ~] = GetPartsNeighbour(I, marksPrev, pointRight(1), pointRight(2), I(pointRight(2),pointRight(1)), ceil(structureTensorRad/3), false);
                    if isempty(marksSelL) || isempty(marksSelR)
                        continue;
                    end
                    
                    ML = mode(marksSelL);
                    MR = mode(marksSelR);
                    
                    % measure distance between ML and MR
                    X_L = LUTable(:, ML);
                    X_R = LUTable(:, MR);
                    dist = sqrt(sum((X_L - X_R).^2));
                    
                    if dist>ignoreThresh
                        curEl = partIDs(ML, MR);
                        outMarks(y,x) = curEl;
                    end
                end
                
 %% Corner parts               
                if cond3  % this is a corner
                    
                    [~, indsX, indsY, ~, ~] = GetPartsNeighbour(I, marksPrev, x, y, I(y,x), ceil(structureTensorRad * 0.6), false);
                    
                    indIm = sub2ind(size(marksPrev), indsY, indsX);
                    depths = I(indIm);
                    points = [indsX', indsY', depths'];
                    pointCentral = [x,y, I(y,x)];
                    pointsLocal = points - repmat(pointCentral, [size(points, 1),1]);
                    partIDs = marksPrev(indIm);
                    
                    % this area should have at least three parts with
                    % partID
                    
                    partIDs = double(partIDs);
                    [elCount, el] = hist(partIDs, unique(partIDs));
                    [elCount, idx] = sort(elCount, 'ascend');
                    el = el(idx);
                    
                    if el > 3
                        el = el(1:3);
                        elCount = elCount(end-2:end);
                    end
                    
                    curr = curr + 1;
                    tablelParts(curr, :) = [el];
                    
                    
                    tableCoords(curr, :) = []; 
                    
                    outMarks(y,x) = curr;
                    
                end
                
            end
        end
        
        %% write the output image

%         imtool(outMarks);
        outMarks = uint16(outMarks);
        imwrite(outMarks, outFile, 'png');
        
        if mod(i,10) == 0
            disp(i);
        end
%         str = ['finish ', num2str(i)];
%         disp(str);

    end


end

function tableParts = recomputePartIds(list_els, lenF, inOut_path, nCurClusters, partIDs, nFirstClusters)
    
    sieveThresh = 0;
    sieveThresh2 = 10;
    
    tableParts = [];
    lenTP = 0;
    partFrequencyAll = zeros(nCurClusters, nCurClusters);
    NumImagesAll = zeros(nCurClusters, nCurClusters);  % in how many images this part expires

    for i = 1:lenF
        
         partFrequencyTemp = zeros(nCurClusters, nCurClusters);
         NumImagesTemp = zeros(nCurClusters, nCurClusters);
         
         I = imread(list_els{i});
         [rows, cols] = find(I>0);
         elsIds = sub2ind(size(I), rows, cols);
         els = I(elsIds);
         els = double(els);
         if length(els) == 0
             continue;
         end
         [elCount, el] = hist(els, unique(els));
         ids = elCount > sieveThresh;
         elCount = elCount(ids);
         el = el(ids);
         
         for j = 1:length(el)
             [lefts, rights] = find(partIDs == el(j));
             left = lefts(1); right = rights(1);
             partFrequencyTemp(left, right) = partFrequencyTemp(left, right) + elCount(j);
             NumImagesTemp(left, right) = NumImagesTemp(left, right) + 1; 
             if left ~= right
                partFrequencyTemp(right, left) = partFrequencyTemp(right, left) + elCount(j);
                NumImagesTemp(right, left) = NumImagesTemp(right, left) + 1;
             end
         end
         partFrequencyAll = partFrequencyAll + partFrequencyTemp;
         NumImagesAll = NumImagesAll + NumImagesTemp;
        if mod(i,200) == 0
            disp(i);
        end
    end
    
    [rows, cols] = find(partFrequencyAll>sieveThresh2);
    Ids = sub2ind(size(partFrequencyAll), rows, cols);
    frequency = partFrequencyAll(Ids); 
    numIm = NumImagesAll(Ids);
    ids = rows<=cols;
    cols = cols(ids);
    rows = rows(ids);
    frequency = frequency(ids);
    numIm = numIm(ids);
    tableParts = [cols, rows, frequency, numIm];
    idsIm = sub2ind(size(partIDs), rows, cols);
    IDs_inImages = partIDs(idsIm);
    [clusterCentres, sampleLables] = ISODATA_3_Layer(tableParts, nFirstClusters);
    % save clusterCentres
    numClusters = size(clusterCentres, 1);
    clusterCentresT = zeros(numClusters, 2);
    for i = 1:numClusters
        clusterCentresT(i, 1) = compute2elementIndex(clusterCentres(i,1), clusterCentres(i,2), nFirstClusters);
        clusterCentresT(i, 2) = compute2elementIndex(clusterCentres(i,3), clusterCentres(i,4), nFirstClusters);
    end
    
    save('clusterTable.mat', 'clusterCentres', 'clusterCentresT', 'numClusters', 'tableParts', 'sampleLables');
    numClusters = size(clusterCentres, 1);
    clusterCentresOut = zeros(numClusters,2);
    for i = 1:numClusters
        clusterCentresOut(i,1) = compute2elementIndex(clusterCentres(i,1), clusterCentres(i,2), nFirstClusters);
        clusterCentresOut(i,1) = compute2elementIndex(clusterCentres(i,3), clusterCentres(i,4), nFirstClusters);
    end
    % create a lookup table
    tableT = zeros(nCurClusters, nCurClusters);
    for i = 1: size(tableParts, 1)
        left = tableParts(i,1); right = tableParts(i,2);
        tableT(left, right) = sampleLables(i);
        tableT(right, left) = sampleLables(i);
    end
    
    % re-compute part ID using this table and WRITE IN PLACE!!!
    parfor i = 1:lenF  
        I = imread(list_els{i});
        [rows, cols] = find(I>0);
        if isempty(rows)
            continue;
        end
        elsIds = sub2ind(size(I), rows, cols);
        el_in_im = I(elsIds);
        elOut = zeros(size(el_in_im));
        for j = 1:length(el_in_im)
            aa = find(IDs_inImages == el_in_im(j));
            if isempty(aa)
                continue;
            end
            aa = aa(1);
            elOut(j) = sampleLables(aa);
        end
        I(elsIds) = elOut;
        disp(list_els{i});
        I = uint8(I);
        imwrite(I, list_els{i}, 'png');
    end
end

function parSave(normFile, Normals)
    save(normFile, 'Normals'); 
end

function statAll = getCooccurrenceStat(list_elsPrev, list_elsCur, listDepth, zScale, structureTensorRad)
    % statistics is collected in the following format: [dx, dy, dz] -
    % offset of the first sub-part w.r.t. the second one

    lenF = length(list_elsPrev);
    load('clusterTable.mat'); % 'clusterCentres', 'clusterCentresT', 'numClusters', 'tableParts', 'sampleLables');
    statAll = {};
    
    for i = 1:size(clusterCentres,1)
        statAll{i} = [];
    end
    
    for i = 1:lenF
        I = imread(listDepth{i});
        I = I(:,:,1);
        elCur = imread(list_elsCur{i});
        elPrev = imread(list_elsPrev{i});
        I = double(I)*zScale;
        
        [rows, cols] = find(elCur > 0);
        if isempty(rows)
            continue;
        end
        ids = sub2ind(size(I), rows, cols);
        els = elCur(ids);
        lenEls = length(rows);
        for j = 1:lenEls
            el = els(j);
            [marksOut, indsXOut, indsYOut, depths, dists] = GetPartsNeighbour(I, elPrev, cols(j), rows(j), I(rows(j), cols(j)), structureTensorRad, false);
            elSP = clusterCentresT(el, :); % subparts
            el_left = min(elSP);
            el_right = max(elSP);
            idsL = find(marksOut ==  el_left);
            idsR = find(marksOut == el_right);
            
            if isempty(idsL) || isempty(idsR)
                continue;
            end
            
            XLeft = indsXOut(idsL); YLeft = indsYOut(idsL); depthLeft = depths(idsL);
            XRight = indsXOut(idsR); YRight = indsYOut(idsR); depthRight = depths(idsR);
            ll = sum([XLeft; YLeft; depthLeft], 2)/length(idsL); 
            rr = sum([XRight; YRight; depthRight] ,2)/length(idsR);
            relPos = rr - ll;
            
            temp = statAll{el};
            temp = [temp, relPos];
            statAll{el} = temp;
        end
        
        disp(i);
    end

end

function plotFrame(V, centre, vecLen, vectColors, n)
    % plot the eigenvectors in 3D
    for ii = 1:n
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [centre(1), centre(1) + vecLen * curVect(1)];
        YY = [centre(2), centre(2) + vecLen * curVect(2)];
        ZZ = [centre(3), centre(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on

    end
end

function plotVector(V, vecLen, vectColors, point)

    XX = [point(1), point(1) + vecLen * V(1)];
    YY = [point(2), point(2) + vecLen * V(2)];
    ZZ = [point(3), point(3) + vecLen * V(3)];

    plot3(XX, YY, ZZ, 'Color', vectColors);
    hold on

end

function visualizeParts(relPos, clusterCentresT, str_folder, cluster1Centres, nClusters, layerID)

    A = exist(str_folder, 'dir');

    if A == 0
        mkdir(str_folder);
    end

    isFIG = false;
    f = figure;
    numParts = size(clusterCentresT, 1);
    fieldSize = [36,36,36];
    fieldCenter = fieldSize / 2;
    elementRadius = 3.5;
    relPos(:,3) = relPos(:,3)/2;
       
    for i = 1:numParts
        disp(i);
        
        positions = [fieldCenter; fieldCenter - relPos(i,:)];
        elements = clusterCentresT(i, :);
        
        surfaceVisualizerT(fieldSize, positions, elements, elementRadius, nClusters, cluster1Centres, []);

        if ~ isFIG
            str = ['layer', num2str(layerID), '_', num2str(i), '.png'];
            str1 = [str_folder, str];
            saveas(f, str1, 'png');
        else
            str = ['layer', num2str(i), '.fig'];
            str1 = [str_folder, str];
            saveas(f, str1, 'fig');
        end
    end
    
    is_ok = true;
end
