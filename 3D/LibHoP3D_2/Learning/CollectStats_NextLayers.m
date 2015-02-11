% this is the main function to co-occurrence statistics for the 4th layer
% and further layers


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatisticsAll, outputCoordsAll, curTSAll] = CollectStats_NextLayers(list_elements, list_depths, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius,...
                                                        nPrevClusters, displacements, lenDisp, layerID, depthStep, dataSetNumber, ...
                                                        is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth, ...
                                                        learningElType, learningElRadius, displ3, cluster1Bounds)
                                                    
    disp('collecting co-occurrence statistics...');
    %  for example displacements = [0 0; 0 -6; 0 6]; or displacements = [0 0; -6, 0; 6, 0];
    %      2 1 3  or
    %     2
    %     1
    %     3
    
    errosionAdder = 0;
    emptyCellID = nPrevClusters + 1;
    isY = true;
    isX = true;
    isX_FB = false;
    
    isTrim = false; % parameters for preliminary processing
    
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end
    
    tic;
    
    outputStatisticsAll = [];
    outputCoordsAll = [];
    curTSAll = 0;
    
    step = 60;
    iiiPrev = 1;
    
    [indsXOut, indsYOut] = getDispAbs(learningElType, learningElRadius);
    
    for iii = 1:step:lenF % This is done to speed up
        
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
                marksPrev = imread(list_elements{i});  % elements of the previous layer

                I = I(:,:,1);

                if dataSetNumber == 2
                    mask = imread(list_mask{i});
                else
                    mask = [];
                end


                if dataSetNumber == 1 || dataSetNumber == 3 % Aim@Shape dataset || Vladislav_STD

                    [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius + errosionAdder, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                            is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                                                                            [],[],[],[]);

                elseif dataSetNumber == 2  % Washington data set


                    [I, Ix, Iy, ~, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius + errosionAdder, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
                                                                            is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                                                                            [],[],[],[]);

                end
                
                Ix = Ix(:,:,1);
                Iy = Iy(:,:,1);
                
                % extend a mask! This is done to treat areas with empty cells
                mask = extendMaskWithDerivatives(mask, cluster1Bounds, Ix, Iy);           

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
                    
                    central = marksPrev(rows(j), cols(j));
                    depthCentral = I(rows(j), cols(j));
                    
                    % check what are left and right neighbours    
                    [indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, cols(j), rows(j), displ3, indsXOut, indsYOut); 

                    % make shure indexes are not out of the image borders
                    [indsXLeft, indsYLeft, indsXRight, indsYRight] = checkImageBoundaries(indsXLeft, indsYLeft, indsXRight, indsYRight, r, c);
                    
                    leftInds  = sub2ind(size(marksPrev), indsYLeft,  indsXLeft);
                    rightInds = sub2ind(size(marksPrev), indsYRight, indsXRight);
                    lefts       = marksPrev(leftInds);
                    rights      = marksPrev(rightInds);
                    leftsMask = mask(leftInds);
                    rightsMask = mask((rightInds));
                    
                    % check is something is empty
                    lenEmpLeft = length(leftsMask(leftsMask == 0))/length(indsYLeft);
                    lenEmpRight = length(rightsMask(rightsMask == 0))/length(indsYRight);
                    
                    if lenEmpLeft >= 0.5 % left should be an empty cell
                        left = emptyCellID;
                        depthLeft = depthCentral;
                    else
                        lefts = lefts(lefts>0);
                        if isempty(lefts)
                            left = 0;
                        else
                            left = mode(lefts);
                            depthLeft =  I(rows(j) + displacements(2,1), cols(j) + displacements(2,2));
                        end

                    end
                    
                    if lenEmpRight >= 0.5 % left should be an empty cell
                        right = emptyCellID;
                        depthRight = depthCentral;
                    else
                        rights = rights(rights>0);
                        if isempty(rights)
                            right = 0;
                        else
                            right = mode(rights);
                            depthRight = I(rows(j) + displacements(3,1), cols(j) + displacements(3,2));   % learn from exact positions
                        end  
                    end



                    if central ~= 0 && left ~= 0 && right ~= 0 && (left ~= emptyCellID || right ~= emptyCellID) % element is detected

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

                        curCoords = [i, cols(j), rows(j)];  % image, x, y
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
    
    disp('Statistics collection time is ...');
    toc

end














