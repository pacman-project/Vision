% this is the main function to co-occurrence statistics for the first and 4th layers
% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatistics, outputCoords, curTS] = CollectStats_3_Layer(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, nClusters, ...
                                                        cluster1Centres, cluster1Lengths, thresh, combs, largestLine, displacements, lenDisp, ... 
                                                        wCoverage, wOverlap, fieldSize, depthStep, is_downsampling, dowsample_rate, dataSetNumber);
    disp('collecting co-occurrence statistics...');
    % field can be non-squared.
    halfFieldSize = floor(fieldSize/2); % for example fieldSize = [17, 5, 71];
    
    if min(halfFieldSize) < 1
        disp('Field size must be larger...');
    end
    
    fieldArea = fieldSize(1) * fieldSize(2);

    shiftX = 0;
    shiftY = 0;
    
    % displacements = [0 0; 0 -6; 0 6];
    %      2 1 3
    
    % should be int8
    outputStatistics = [];
    outputCoords = [];
    curTS = 0;
    
    isY = true;
    isX = true;
    isTrim = false; % parameters for preliminary processing

    parfor i = 1:lenF 
        
        I = imread(list_depth{i});
        I = I(:,:,1);
        mask = [];
        
        if dataSetNumber == 1  % Aim@Shape dataset
            if is_downsampling
                I = imresize(I, dowsample_rate);
            end
            
            [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discSize, isX, isY, ...
                                isTrim, dxKernel, sigmaKernelSize, sigma);
 
                            
        elseif dataSetNumber == 2  % Washington data set
            
            mask = imread(list_mask{i});
            if is_downsampling
                I = imresize(I, dowsample_rate);
                mask = imresize(mask, dowsample_rate);
            end           
            [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, isX, isY,...
                                isTrim, dxKernel, sigmaKernelSize, sigma);
        end
        
        if ~is_successfull
            continue;
        end
        
        % fieldsize and center are ordered:      x, y and z directions
        
        windowMarks = zeros(fieldSize(2), fieldSize(1));
        windowErrors1 = zeros(fieldSize(2), fieldSize(1)); % here we collect errors for each mark
        windowMarksAlt = zeros(fieldSize(2), fieldSize(1));
        windowErrors1Alt = zeros(fieldSize(2), fieldSize(1)); % altermative marks and errors


        for j = halfFieldSize(2)+1 : 6 : r-halfFieldSize(2)     % y-direction (rows)
            for k = halfFieldSize(1)+1 : 6 : c-halfFieldSize(1)   % x-direction (columns)
                
                % working only with locations where the full element can be
                % detected
                windowMask = mask(j-halfFieldSize(2) : j+halfFieldSize(2), k-halfFieldSize(1) : k+halfFieldSize(1));
                
                if length(windowMask(windowMask == 1)) < 0.9 * fieldArea % no mask in the field of view
                    continue;
                end
                
                windowI =  I(j-halfFieldSize(2)  : j+halfFieldSize(2), k-halfFieldSize(1) : k+halfFieldSize(1));
                windowIx = Ix(j-halfFieldSize(2) : j+halfFieldSize(2), k-halfFieldSize(1) : k+halfFieldSize(1));
                windowIy = Iy(j-halfFieldSize(2) : j+halfFieldSize(2), k-halfFieldSize(1) : k+halfFieldSize(1));

                
                % discretize every other row
                for jj = 1:2:fieldSize(2)
                    
                    inds = 3:6:fieldSize(1)-1;
                    tempLine = windowIx(jj,:);
                    fx = tempLine(inds);
                    strLen = length(fx);
                    [nearestClusters, errors, alternativeClusters, alternativeErrors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, cluster1Lengths, thresh);
                    % [output, curErrs] = lineDiscretizationOptLayer1(nearestClusters, errors,  wCoverage, wOverlap, combs, largestLine);
                    
                    tempLine = zeros(1, fieldSize(1));
                    tempLine(inds) = nearestClusters;
                    windowMarks(jj,:) = tempLine;
                    
                    tempLine = zeros(1, fieldSize(1));
                    tempLine(inds) = errors;
                    windowErrors1(jj,:) = tempLine;
                    
                    tempLine = zeros(1, fieldSize(1));
                    tempLine(inds) = alternativeClusters;
                    windowMarksAlt(jj,:) = tempLine;
                    
                    tempLine = zeros(1, fieldSize(1));
                    tempLine(inds) = alternativeErrors;
                    windowErrors1Alt(jj,:) = tempLine;

                end
                
%                 windowMarks =   windowMarks.*windowMask; 
%                 windowErrors1 = windowErrors1.*windowMask;
%                 windowMarksAlt = windowMarksAlt.*windowMask; 
%                 windowErrors1Alt = windowErrors1Alt.*windowMask;

                % next we create second layer elements in every window
                
                center = halfFieldSize + [1,1,1];         

                isActive = zeros(1,lenDisp);
                line = zeros(1, lenDisp*2 - 1); 
                line = int8(line);
                curCoords = zeros(1,3);
                cont = true;
                
                % displacements are ordered like that:  (y,x)
                % fieldsize and center are ordered:      x, y and z directions

                for jj = 1:lenDisp
                    if cont  % this is done instead of writing break loop (which is not allowed by parfor)
                        
                        % here we extract lines of size 3 in y direction
                        curDisp = displacements(jj, :);
                        curCenter = [center(2) + curDisp(1), center(1) + curDisp(2)];
                        sliceMarks1 = windowMarks(curCenter(1) - halfFieldSize(2): curCenter(1) + halfFieldSize(2), curCenter(2):curCenter(2));
                        sliceErrors1 = windowErrors1(curCenter(1) - halfFieldSize(2): curCenter(1) + halfFieldSize(2), curCenter(2):curCenter(2));
                        sliceMarks1Alt = windowMarksAlt(curCenter(1) - halfFieldSize(2): curCenter(1) + halfFieldSize(2), curCenter(2):curCenter(2));
                        sliceErrors1Alt = windowErrors1Alt(curCenter(1) - halfFieldSize(2): curCenter(1) + halfFieldSize(2), curCenter(2):curCenter(2));
                        
%                         indsY = 1:2:fieldSize(2);
%                         sliceMarks1 = sliceMarks1(indsY);
%                         sliceErrors1 = sliceErrors1(indsY);
%                         sliceMarks1Alt = sliceMarks1Alt(indsY);
%                         sliceErrors1Alt = sliceErrors1Alt(indsY);
                        
                        [nearestCluster, error] = discretizeY(sliceMarks1, sliceErrors1, sliceMarks1Alt, sliceErrors1Alt, 5, 0);
                        %[nearestCluster, error] = discretizeSlice(sliceMarks1, sliceErrors1, 3);

                        if nearestCluster(3) > 0 && error(3) < 1
                            % define 2nd layer element
                            dY = windowIy(curCenter(1), curCenter(2));
                            clusterY = define1Cluster(dY, nClusters, cluster1Lengths, thresh);
                            curEl = compute2elementIndex(nearestCluster(3), clusterY, nClusters);
                            isActive(jj) =  1; % redundant but just to check

                            % now write this to the array line
                            if jj == 1
                                depthCentral = windowI(curCenter(1), curCenter(2));
                                line(1) = curEl;
                            else
                                curDepth = windowI(curCenter(1), curCenter(2));
                                relDepth = (curDepth - depthCentral) / depthStep; % relative depth
                                if relDepth > 127
                                    relDepth = 127;
                                elseif relDepth < -127
                                    relDepth = -127; % for one bite coding
                                end
                                line((jj-1)*2) = curEl;
                                line((jj-1)*2+1) = relDepth;
                            end
                        else
                            isActive(jj) =  0;
                            cont = false;
                        end
                    end
                end
                
                if sum(isActive) == lenDisp
                    curCoords = [i, k, j];
                    curCoords = uint16(curCoords);
                    outputStatistics = [outputStatistics; line];
                    outputCoords = [outputCoords; curCoords];
                    curTS = curTS + 1;
                end

            end
        end

    i
    end

end












