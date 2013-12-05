% this is the main function to co-occurrence statistics for the first and 4th layers
% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatistics, curTS] = CollectStats_3_4_Layers(list_depth, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, cluster1Centres, cluster1Lengths, thresh, ...
                 outRoot, combs, largestLine, displacements, wCoverage, wOverlap, fieldSize, depthStep)

    if fieldSize(1) == fieldSize(2)
        halfFieldSize = floor(fieldSize(1)/2);
    else
        halfFieldSize = floor(fieldSize/2);
        disp('Something is wrong...');
    end
    fieldArea = fieldSize(1) * fieldSize(2);

    shiftX = 0;
    shiftY = 0;
    
    %        2       
    %      5 1 3      
    %        4
    % displacements = [0,0; 0,-4; 4,0; 0,4; -4,0]; 

    %      9 2 6      
    %      5 1 3      
    %      8 4 7
    % displacements = [0,0; 0,-4; 4,0; 0,4; -4,0; 4,-4; 4,4; -4,4; -4,-4];  % should work
    lenDisp = size(displacements, 1);
    
    % should be int8
    outputStatistics = [];
    curTS = 0;
    
    isErrosion = false;
    discSize = 0;
    isY = true;
    isTrim = true;

    parfor i = 1:lenF % To use later for models
        
        I = imread(list_depth{i});

        %this function returns the trimmed image and the mask
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, isErrosion, discSize, isY, isTrim, dxKernel, sigmaKernelSize, sigma);

        if ~is_successfull
            continue;
        end
        
        windowMarks = zeros(fieldSize(1), fieldSize(2));
        windowErrors1 = zeros(fieldSize(1), fieldSize(2)); % here we collect errors for each mark
        windowI = zeros(fieldSize(1), fieldSize(2));
        windowIx = zeros(fieldSize(1), fieldSize(2));
        windowIy = zeros(fieldSize(1), fieldSize(2));
        windowMask = zeros(fieldSize(1), fieldSize(2));


        for j = halfFieldSize+1 : 4 : r-halfFieldSize 
            % we now put every element to the central position of the windows
            for k = halfFieldSize+1 : 4 : c-halfFieldSize
                
                % working only with locations where the full element can be
                % detected
                windowMask = mask(j-halfFieldSize : j+halfFieldSize, k-halfFieldSize : k+halfFieldSize);
                
                if length(windowMask(windowMask == 1)) < fieldArea % no mask in the field of view
                    continue;
                end
                
                windowI =  I(j-halfFieldSize : j+halfFieldSize, k-halfFieldSize : k+halfFieldSize);
                windowIx = Ix(j-halfFieldSize : j+halfFieldSize, k-halfFieldSize : k+halfFieldSize);
                windowIy = Iy(j-halfFieldSize : j+halfFieldSize, k-halfFieldSize : k+halfFieldSize);
                
                % discretize every row
                for jj = 1:fieldSize(1)
                    inds = 3:2:fieldSize(1)-1;
                    line = windowIx(jj,:);
                    fx = line(inds);
                    strLen = length(fx);
                    [nearestClusters, errors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, cluster1Lengths, thresh);
                    [output, curErrs] = lineDiscretizationOptLayer1(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine);
                    line = zeros(1, fieldSize(1));
                    line(inds) = output;
                    windowMarks(jj,:) = line;
                    
                    line = zeros(1, fieldSize(1));
                    line(inds) = curErrs;
                    windowErrors1(jj,:) = line;
                end
                windowMarks =   windowMarks.*windowMask; 
                windowErrors1 = windowErrors1.*windowMask; 

                % next we create second layer elements in every window
                % we predict centers in the positions:  3,7,11 (if fieldsize == 13)
                
                center = halfFieldSize+1;
                % displacements = [0,0; 0,-4; 4,0; 0,4; -4,0; 4,-4; 4,4; -4,4; -4,-4];;
                isActive = zeros(1,lenDisp);
                line = zeros(1, lenDisp*2 - 1); 
                line = int8(line);
                cont = true;

                for jj = 1:lenDisp
                    if cont  % this is done instead of writing break loop (which is not allowed by parfor)
                        curDisp = displacements(jj, :);
                        curCenter = [center + curDisp(1), center + curDisp(2)];
                        sliceMarks1 = windowMarks(curCenter(1) - 1: curCenter(1) + 1, curCenter(2) - 2:curCenter(2) + 2);
                        sliceErrors1 = windowErrors1(curCenter(1) - 1: curCenter(1) + 1, curCenter(2) - 2:curCenter(2) + 2);
                        [nearestCluster, error] = discretizeSlice(sliceMarks1, sliceErrors1, 3);

                        if nearestCluster(2) > 0 && error(2) < 1
                            % define 2nd layer element
                            dY = windowIy(curCenter(1), curCenter(2));
                            clusterY = define1Cluster(dY, nClusters, cluster1Lengths, thresh);
                            curEl = compute2elementIndex(nearestCluster(2), clusterY, nClusters);
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
                
                if sum(isActive) == 9
                    outputStatistics = [outputStatistics; line];
                    curTS = curTS + 1;
                end

            end
        end

    i
    end

end












