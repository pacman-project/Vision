% this is the main function to co-occurrence statistics for the 4th layer
% and further layers


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatisticsAll, outputCoordsAll, curTSAll] = CollectStats_NextLayers(list_elements, list_depths, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                                                        displacements, lenDisp, displacement34, layerID, depthStep, dataSetNumber, ...
                                                        is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth)
                                                    
    disp('collecting co-occurrence statistics...');
    %  for example displacements = [0 0; 0 -6; 0 6]; or displacements = [0 0; -6, 0; 6, 0];
    %      2 1 3  or
    %     2
    %     1
    %     3
    
    
    isY = false;
    isX = false;
    isTrim = false; % parameters for preliminary processing
    
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end
    
    tic;
    
    outputStatisticsAll = [];
    outputCoordsAll = [];
    curTSAll = 0;
    
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
        
    
        parfor i = indStart:indEnd

                I = imread(list_depths{i});     % depth image
                marks3 = imread(list_elements{i});  % elements of the previous layer

                I = I(:,:,1);

                if dataSetNumber == 2
                    mask = imread(list_mask{i});
                else
                    mask = [];
                end


                if dataSetNumber == 1 || dataSetNumber == 3 % Aim@Shape dataset || Vladislav_STD

                    [I, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                            is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);

                elseif dataSetNumber == 2  % Washington data set


                    [I, ~, ~, ~, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                            is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);

                end

                if ~is_successfull
                    continue;
                end

                [rEl, cEl] = size(mask);  % all three images should be of the same size!
                if r~=rEl || c ~= cEl
                    disp('ERROR');
                end
                [rEl, cEl] = size(marks3);
                if r~=rEl || c ~= cEl
                    disp('ERROR');
                end
                if ~is_successfull
                    disp('ERROR');
                end

                % next we try to match the next layer element around each window


                [rows, cols] = find(marks3 > 0);
                nEl = length(rows);

                for j = 1: nEl

                    if layerID == 4 || layerID == 6 || layerID == 8
                        if rows(j) < displacement34 + 1  || rows(j) > r - displacement34  || cols(j) < displacement34 + 1 || cols(j) > c - displacement34
                            continue;
                        end
                    elseif layerID == 3 || layerID == 5 || layerID == 7 
                        if cols(j) < displacement34 + 1 || cols(j) > c - displacement34
                            continue;
                        end
                    end

                    central = marks3(rows(j) + displacements(1,1), cols(j) + displacements(1,2));
                    depthCentral = I(rows(j) + displacements(1,1), cols(j) + displacements(1,2));
                    left =  marks3(rows(j) + displacements(2,1), cols(j) + displacements(2,2)); 
                    right = marks3(rows(j) + displacements(3,1), cols(j) + displacements(3,2));  
                    depthLeft =  I(rows(j) + displacements(2,1), cols(j) + displacements(2,2));
                    depthRight = I(rows(j) + displacements(3,1), cols(j) + displacements(3,2));   % learn from exact positions


                    if central ~= 0 && left ~= 0 && right ~= 0  % element is detected

                        line = zeros(1, lenDisp*2 - 1); 
                        line(1) = central;
                        line(2) = left;
                        line(4) = right;

                        relDepthL = round((depthLeft - depthCentral) / depthStep);  % relative depth
                        if relDepthL > maxRelDepth
                            relDepthL = maxRelDepth;
                        elseif relDepthL < -maxRelDepth
                            relDepthL = -maxRelDepth; % for one bite coding
                        end
                        line(3) = relDepthL;

                        relDepthR = round((depthRight - depthCentral) / depthStep);  % relative depth
                        if relDepthR > maxRelDepth
                            relDepthR = maxRelDepth;
                        elseif relDepthR < -maxRelDepth
                            relDepthR = -maxRelDepth; % for one bite coding
                        end
                        line(5) = relDepthR;

                        line = int16(line);

                        curCoords = [i, cols(j) + displacements(1,2), rows(j) + displacements(1,1)];  % image, x, y
                        curCoords = uint16(curCoords);
                        outputCoords = [outputCoords; curCoords];
                        outputStatistics = [outputStatistics; line];
                        curTS = curTS + 1;

                    end
                end

            i
        end
        
        outputStatisticsAll = [outputStatisticsAll; outputStatistics];
        outputCoordsAll = [outputCoordsAll; outputCoords];
        curTSAll = curTSAll + curTS;

        iiiPrev = iii + 1;
    
    end
    
    toc

end














