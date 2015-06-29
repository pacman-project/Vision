% this is the main function to co-occurrence statistics for the 4th layer
% and further layers


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatisticsAll, outputCoordsAll, outputScalesAll, outputFramesAll, curTSAll] = CollectStats_NextLayers(list_els, list_input, list_mask, lenF, filtOptions,  ...
                                                        nPrevClusters, layerID, depthStep, dataSetNumber, inputDataType, zScale, ...
                                                        maxRelDepth, cluster1Bounds, nClusters, options, offsetsConventional)

% (list_elements, list_depths, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius,...
%                                                         nPrevClusters, displacements, lenDisp, layerID, depthStep, dataSetNumber, ...
%                                                         is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth, ...
%                                                         learningElType, learningElRadius, displ3, cluster1Bounds)
                                                    
    is_visualization = false;                                                
                                                    
    disp('collecting co-occurrence statistics...');
    %  for example displacements = [0 0; 0 -6; 0 6]; or displacements = [0 0; -6, 0; 6, 0];
    %      2 1 3  or
    %     2
    %     1
    %     3

    emptyCellID = nPrevClusters + 1;
    
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end
    
    tic;
    
    outputStatisticsAll = [];
    outputCoordsAll = [];
    outputScalesAll = [];
    outputFramesAll = [];
    curTSAll = 0;
    
    step = 60;
    iiiPrev = 1;
    
    maxAngle = options.maxAngleToZaxis;
    
%     [indsXOut, indsYOut] = getDispAbs(learningElType, learningElRadius);
    
    for iii = 1:step:lenF % This is done to speed up
        
        outputStatistics = [];
        outputCoords = [];
        outputScales = [];
        outputFrames = [];
        curTS = 0;
        
        indStart = iiiPrev;
        rest = lenF - iii;
        if rest < step
            indEnd = lenF;
        else
            indEnd = iii;
        end
    
        parfor i = indStart:indEnd
                      
%                 
%                 if ~strcmp(list_input{i},'D:/Input Data/VladislavSTD/Vladislav_STD/depth/99_3_2_2_3_2.png')
%                     continue;
%                 end
                numEl = 0;
                
%                 i
                
                I = imread(list_input{i});     % depth image
                marksPrev = imread(list_els{i});  % elements of the previous layer
                % read image with normals
                normFile = [list_els{i}(1:end - 4), '_N.mat']; 
                Normals = load(normFile); 
                Normals = Normals.Normals;

                I = I(:,:,1);
                I = I * zScale;
                [r,c] = size(I);

                if dataSetNumber == 2
                    mask = imread(list_mask{i});
                else
                    mask = [];
                end

                if dataSetNumber == 1 || dataSetNumber == 3 || dataSetNumber == 4 % Aim@Shape dataset || Vladislav_STD

%                     [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius + errosionAdder, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
%                                                                             is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, ...
%                                                                             [],[],[],[]);

                    [I, ~, ~, mask, is_successfull] = preliminaryProcessing(I, mask, filtOptions);                                                    

                elseif dataSetNumber == 2  % Washington data set

% 
%                     [I, Ix, Iy, ~, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius + errosionAdder, isX, isY, isX_FB, isTrim, dxKernel, sigmaKernelSize, sigma, ...
%                                                                             is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, ...
%                                                                             [],[],[],[]);
                    [I, ~, ~, mask, is_successfull] = preliminaryProcessing(I, mask, filtOptions);

                end
                
%                 Ix = Ix(:,:,1);
%                 Iy = Iy(:,:,1);
%                 
%                 % extend a mask! This is done to treat areas with empty cells
%                 mask = extendMaskWithDerivatives(mask, cluster1Bounds, Ix, Iy);           

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

                % next we try to match the next layer element around each
                % window
                
                marksPrev(mask == 0) = 0;
                
                [rows, cols] = find(marksPrev > 0);
                nEl = length(rows);
                
                

                for j = 1: nEl 
                    
%                     disp(j);
%                     if j ~= 138
%                         continue;
%                     end
                                    
                    y = rows(j); 
                    x = cols(j);
                    
%                     if x~= 305 || y~= 127
%                         continue;
%                     end
                    
                    depthCentral = I(y, x);
                    Norm = Normals(y,x,:);
                    Norm = squeeze(Norm);
                    zAxis = [0;0;1];
                    angle = acos(dot(Norm, zAxis))*180/pi;
                    
                                        
                    if angle > maxAngle % The slope is too high to estimate robustly!
                        continue;
                    end
                    
                    if layerID == 3
                        central = compute2elementIndex(ceil(nClusters/2), ceil(nClusters/2), nClusters);
                        elScale = marksPrev(y,x);
                        is2D = false;
                    end
                    
                    if layerID == 3 || layerID == 4
                        
                        receptiveFieldRad = elScale * 4 + 1;
                        receptiveFieldRad = double(receptiveFieldRad);
                        
                        % estimate a darboux frame at this point   
                        [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, receptiveFieldRad, is2D, filtOptions.minDepth, filtOptions.maxDepth);

                        if length(indsXOut) < options.ignoreThresh{layerID}
                            continue;
                        end

                        % now we convert these coordinates to the local frame of
                        % reference (centered at the point [j,i,curDepth])
                        xs = indsXOut - x;
                        ys = indsYOut - y;
                        zs = depths - I(y,x);

                        if options.methodId == 1  % paraboloid fitting                  
                            [H, K, QF1, QF2, S, V, D, is_ok] = paraboloidFitting(xs, ys, zs);
                        elseif options.methodId == 2 % method of Gabriel Taubin (modified)
                            disp('Sorry Gabriel Taubins method is not implemented');
                        end
                        
                        if ~is_ok
                            continue;
                        end
                        
                        if layerID == 3
                            curDisp = V(:,2);
                        elseif layerID == 4
                            curDisp = V(:,3);
                        end
                    end
                    
                    offset = receptiveFieldRad;
                    curOffset = zeros(3, 2);
                    curOffset(:, 1) =   curDisp * offsetsConventional{layerID}; 
                    curOffset(:, 2) = - curDisp * offsetsConventional{layerID};
                    
                    % find a surface point that is closest to the [x,y,depthCentral]
                    [xLeft, yLeft, zLeft] = findClosestPoint(I, [x + curOffset(1, 1), y + curOffset(2, 1),depthCentral + curOffset(3, 1)]);
                    [xRight, yRight, zRight] = findClosestPoint(I, [x + curOffset(1, 2), y + curOffset(2, 2),depthCentral + curOffset(3, 2)]);
                    
                    if is_visualization
                        shift = 30;
                        centre = shift+1;
                        II = I(y-30: y+30, x - 30: x+30);
                        II(II < 1050) = NaN;
                        % visualize this image in 3D
                        surf(II, 'FaceColor',[0.3, 0.3 0.3], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceLighting', 'phong');
                        camlight left
                        axis equal;
                        hold on
                        curDepth = II(shift+1,shift+1);
                        
                        vectColors = eye(3);
                        vecLen = 20;
                        plotFrame(V, [centre, centre, depthCentral], vecLen, vectColors, 3);
                        
                        hold on
                        plotFrame(eye(3), [centre, centre, depthCentral], vecLen, vectColors, 3);
                        
                        hold on
                        XXX = xs + centre;
                        YYY = ys + centre;
                        ZZZ = zs + depthCentral;
                        scatter3(XXX, YYY, ZZZ, 'red');
                        hold on

                    end
                    

                    % find those parts that are detected at these offsets
                    [lefts, indsXLeft, indsYLeft, depthsLeft, distsLeft] = GetPartsNeighbour(I, marksPrev, xLeft, yLeft, zLeft, floor(offsetsConventional{layerID}/2), is2D);
                    [rights, indsXRight, indsYRight, depthsRight, distsRight] = GetPartsNeighbour(I, marksPrev, xRight, yRight, zRight, floor(offsetsConventional{layerID}/2), is2D);
                    
                    if isempty(lefts) || isempty(rights)
                        continue;
                    end
 
                    if layerID == 3
                        
                        if length(lefts) > 1  % find the one with the largest scale and shortest distance                       
                           lenE = length(lefts);
                           score = double(options.weight * lefts) - distsLeft;
                           id = find(score == max(score));
                           if length(id) > 1
                               id = id(1);
                           end
                           [lefts, indsXLeft, indsYLeft, distsLeft] = myFilter4(lefts,indsXLeft,indsYLeft,distsLeft, id);
                           depthsLeft= depthsLeft(id);
                        end

                        if length(rights) > 1  % find the one with the largest scale and shortest distance 
                           lenE = length(rights);
                           score = double(options.weight * rights) - distsRight;
                           id = find(score == max(score));
                           if length(id) > 1
                               id = id(1);
                           end
                           [rights, indsXRight, indsYRight, distsRight] = myFilter4(rights,indsXRight,indsYRight,distsRight, id);
                           depthsRight = depthsRight(id);
                        end
                        
                        if is_visualization
                            XXX = [indsXLeft, indsXRight] - x + centre;
                            YYY = [indsYLeft, indsYRight] - y + centre;
                            ZZZ = [depthsLeft; depthsRight];
                            scatter3(XXX, YYY, ZZZ, 'green');
                            hold on
                            
                        end
                        
                        curScale = min([lefts, rights, elScale]);
                        
                        NormalLeft = squeeze(Normals(indsYLeft, indsXLeft, :));
                        NormalRight = squeeze(Normals(indsYRight, indsXRight, :));
                        NormalCentral = Norm;
                        
                        if is_visualization
%                             plotFrame(NormalLeft, [indsXLeft-x+centre, indsYLeft-y+centre, depthsLeft], vecLen, vectColors, 1);
%                             plotFrame(NormalRight,[indsXRight-x+centre, indsYRight-y+centre, depthsRight], vecLen, vectColors, 1);
                        end
                        
                        % compute and quantize relative orienations of
                        % normals
                        Xtemp = V(:,2);
                        Ytemp = V(:,3);  % local x and y axis
                        
                        % compute projections
                        angleXLeft =  90 - 180 * acos(dot(Xtemp,NormalLeft)) /pi;
                        angleYLeft =  90 - 180 * acos(dot(Ytemp,NormalLeft)) /pi;
                        angleXRight = 90 -  180 * acos(dot(Xtemp,NormalRight))/pi;
                        angleYRight = 90 - 180 * acos(dot(Ytemp,NormalRight))/pi;
                        
                        clusterXLeft = define1Cluster(angleXLeft, cluster1Bounds, nClusters);
                        clusterYLeft  = define1Cluster(angleYLeft, cluster1Bounds, nClusters);
                        clusterXRight = define1Cluster(angleXRight, cluster1Bounds, nClusters);
                        clusterYRight = define1Cluster(angleYRight, cluster1Bounds, nClusters);
                        
                        if clusterXLeft <= 0 || clusterYLeft <= 0 || clusterXRight <=0 || clusterYRight <=0
                            continue;
                        end
                        
                        % define left and right neighbours
                        left = compute2elementIndex(clusterXLeft, clusterYLeft, nClusters);
                        right = compute2elementIndex(clusterXRight, clusterYRight, nClusters);
                        a = 1;
                        
                        if left == 47 && right == 40
                            a = 2;
                        elseif left == 34 && right == 20
                            a = 2;
                        elseif left == 30 && right == 16
                            a = 2;
                        elseif left == 23 && right == 23
                            a = 2;
                        end
                           
                        if a == 2
                            a = 3;
                        end
                        
                        % compute relative depths of the pixels
                        T = eye(4); T(1,4) = -x; T(2,4) = -y; T(3,4) = -depthCentral;
                        R = eye(4,4); R(1:3, 1) = Xtemp'; R(1:3, 2) = Ytemp'; R(1:3, 3) = NormalCentral';
                        M = inv(R)*T; 
                        
                        % coordinates of the neighbours in the local frame
                        % of reference
                        coordsLeft = M * [indsXLeft;  indsYLeft;  depthsLeft;  1];
                        coordsRight = M * [indsXRight; indsYRight; depthsRight; 1];
                        
                        offsetXLeft = abs(coordsLeft(1));
                        offsetXRight = abs(coordsRight(1));
                        
                        depthLeft = coordsLeft(3)  * offsetsConventional{layerID}/offsetXLeft;
                        depthRight = coordsRight(3)* offsetsConventional{layerID}/offsetXRight;
                        
                        if offsetXLeft == 0 || offsetXRight == 0;
                            continue;
                        end
                    end
                    


                    if central ~= 0 && left ~= 0 && right ~= 0 && (left ~= emptyCellID || right ~= emptyCellID) % element is detected

                        line = zeros(1, 3*2 - 1); 
                        line(1) = central;
                        line(2) = left;
                        line(4) = right;

                        relDepthL = round(depthLeft / depthStep);  % relative depth
                        if relDepthL > maxRelDepth
                            relDepthL = maxRelDepth;
                        elseif relDepthL < -maxRelDepth
                            relDepthL = -maxRelDepth; % for one bite coding
                        end
                        line(3) = relDepthL;

                        relDepthR = round(depthRight / depthStep);  % relative depth
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
                        scaleLine = [offsetXLeft, offsetXRight, curScale];
                        outputScales = [outputScales; scaleLine];
                        outputFrames = [outputFrames; V(:)'];
                        
                        curTS = curTS + 1;
                        numEl = numEl +1;

                    end
                end
%             if numEl == 0;
%                 list_input{i}
%             end
            i
        end
        
        outputStatisticsAll = [outputStatisticsAll; outputStatistics];
        outputCoordsAll = [outputCoordsAll; outputCoords];
        outputScalesAll = [outputScalesAll; outputScales];
        outputFramesAll = [outputFramesAll; outputFrames];
        curTSAll = curTSAll + curTS;

        iiiPrev = iii + 1;
    
    end
    
    disp('Statistics collection time is ...');
    toc

end

function [af, bf, cf, df] = myFilter4(a,b,c,d, ids)
    af = a(ids);
    bf = b(ids);
    cf = c(ids);
    df = d(ids);
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





%                     leftInds  = sub2ind(size(marksPrev), indsYLeft,  indsXLeft);
%                     rightInds = sub2ind(size(marksPrev), indsYRight, indsXRight);
%                     lefts       = marksPrev(leftInds);
%                     rights      = marksPrev(rightInds);
%                     leftsMask = mask(leftInds);
%                     rightsMask = mask((rightInds));
%                     
%                     % check is something is empty
%                     lenEmpLeft = length(leftsMask(leftsMask == 0))/length(indsYLeft);
%                     lenEmpRight = length(rightsMask(rightsMask == 0))/length(indsYRight);
%                     
%                     if lenEmpLeft >= 0.5 % left should be an empty cell
%                         left = emptyCellID;
%                         depthLeft = depthCentral;
%                     else
%                         lefts = lefts(lefts>0);
%                         if isempty(lefts)
%                             left = 0;
%                         else
%                             left = mode(lefts);
%                             coordsLeft = [min(rows(j) + displacements(2,1), r), min(cols(j) + displacements(2,2),c)];
%                             depthLeft =  I(coordsLeft(1), coordsLeft(2));
%                         end
% 
%                     end
%                     
%                     if lenEmpRight >= 0.5 % left should be an empty cell
%                         right = emptyCellID;
%                         depthRight = depthCentral;
%                     else
%                         rights = rights(rights>0);
%                         if isempty(rights)
%                             right = 0;
%                         else
%                             right = mode(rights);
%                             coordRight = [min(rows(j) + displacements(3,1), r), min(cols(j) + displacements(3,2), c)];
%                             depthRight = I(coordRight(1), coordRight(2));   % learn from exact positions
%                         end  
%                     end













