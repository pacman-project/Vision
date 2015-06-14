% this function collects co-occurrence statistics each layer starting form
% the 3rd one and encode it as 3D statistical maps


function [statistics, sumSamples] = buildStatMap_NextLayers(list_elements, list_depths, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius,...
                                                        nPrevClusters, depthStep, dataSetNumber, ...
                                                        is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, fieldSize, filterThresh)
                                                    
    disp('constructing statistical maps...');
    %  for example displacements = [0 0; 0 -6; 0 6]; or displacements = [0 0; -6, 0; 6, 0];
    %      2 1 3  or
    %     2
    %     1
    %     3
    
    errosionAdder = 0;
    emptyCellID = nPrevClusters + 1;
    isY = false;
    isX = false;
    isX_FB = false;
    
    halfFieldSize = floor(fieldSize/2);
    maxRelDepth = halfFieldSize(3);
    
    
    isTrim = false; % parameters for preliminary processing
    
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end
    
    tic;
    statBox = zeros(fieldSize(1), fieldSize(2), fieldSize(3));
    
    for i = 1:nPrevClusters
        for j = 1:nPrevClusters
            statistics{i}{j} = statBox;
        end
    end
    
    
    sumSamples = zeros(nPrevClusters, nPrevClusters);
    halfFieldSize = floor(fieldSize/2);
    fieldCentre = ceil(fieldSize/2);
    
    
    step = 100;
    iiiPrev = 1;
    
    for iii = step:step:lenF % This is done to speed up
        
        outputStatistics = [];
        outputCoords = [];
        curTS = 0;
        
        indStart = iiiPrev;
        rest = lenF - iii;
        if rest < step
            indEnd = lenF;
        else
            indEnd = iii;
        end
        
    
        for i = indStart:indEnd

                I = imread(list_depths{i});     % depth image
                marksPrev = imread(list_elements{i});  % elements of the previous layer

                I = I(:,:,1);

                if dataSetNumber == 2
                    mask = imread(list_mask{i});
                else
                    mask = [];
                end


                if dataSetNumber == 1 || dataSetNumber == 3 % Aim@Shape dataset || Vladislav_STD

                    [I, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius + errosionAdder, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                            is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                                                                            [],[],[],[]);

                elseif dataSetNumber == 2  % Washington data set


                    [I, ~, ~, ~, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius + errosionAdder, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                            is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                                                                            [],[],[],[]);

                end

                if ~is_successfull
                    continue;
                end

                [rEl, cEl] = size(mask);  % all three images should be of the same size!
                if r~=rEl || c ~= cEl
                    disp('ERROR');
                end
                [rEl, cEl] = size(marksPrev);
                if r~=rEl || c ~= cEl
                    disp('ERROR');
                end
                if ~is_successfull
                    disp('ERROR');
                end

                % next we try to match the next layer element around each window


                [rows, cols] = find(marksPrev > 0);
                nEl = length(rows);

                for j = 1: nEl


                    el = marksPrev(rows(j), cols(j));
                    depthCentral = I(rows(j), cols(j));
                    
                    indsW = find(rows >= rows(j) - halfFieldSize(2) & rows <= rows(j) + halfFieldSize(2) & cols >= cols(j) - halfFieldSize(1) & cols <= cols(j) + halfFieldSize(1));
                    
                    lenIW = length(indsW);
                    
                    for iiii = 1:lenIW
                        
                        curEl = marksPrev(rows(indsW(iiii)), cols(indsW(iiii)));
                        relDepth = round((I(rows(indsW(iiii)), cols(indsW(iiii))) - depthCentral) / depthStep);  % relative depth
                        shiftX = fieldCentre(1) + cols(indsW(iiii)) - cols(j);
                        shiftY = fieldCentre(2) + rows(indsW(iiii)) - rows(j);
                        
                        if relDepth > maxRelDepth
                            relDepth = maxRelDepth;
                        elseif relDepth < -maxRelDepth
                            relDepth = -maxRelDepth; % for one bite coding
                        end
                        relDepth = relDepth + fieldCentre(3);
                        
                        
                        
                        
                        
                        statistics{el}{curEl}(shiftX, shiftY, relDepth) = statistics{el}{curEl}(shiftX, shiftY, relDepth) + 1;
                        sumSamples(el, curEl) = sumSamples(el, curEl) + 1;
                    end


                end

            i
        end

        iiiPrev = iii + 1;
    
    end
    
    % filter out those pairs of statistics which 
    
    for i = 1:nPrevClusters
        for j = 1:nPrevClusters
            if sumSamples(i, j) < filterThresh
                statistics{i}{j} = [];
            end
        end
    end

end














