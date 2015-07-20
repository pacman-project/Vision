% this is to learn co-occurrence statistics from meshes and point clouds


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatisticsAll, outputCoordsAll, outputFramesAll, curTSAll] = CollectStats_NextLayersMesh_PC(list_els, list_input, lenF, ...
                                                        nPrevClusters, layerID, depthStep, dataSetNumber, ...
                                                        maxRelDepth, cluster1Bounds, nClusters, offsetsConventional)

                                                    
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
    nf = 20;
    
    outputStatisticsAll = [];
    outputCoordsAll = [];
    outputScalesAll = [];
    outputFramesAll = [];
    curTSAll = 0;
    
    step = 60;
    iiiPrev = 1;
    
    if layerID == 3;
        central = compute2elementIndex(ceil(nClusters/2), ceil(nClusters/2), nClusters);
    end
    
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
    
        for i = 2:2%indStart:indEnd
                      
            numEl = 0;

            fileName = list_input{i};
            outFileM = [fileName(1:end-4), '.mat'];   
            outFileAP = [fileName(1:end-4), 'AP.mat']; 
            
            [V, F, ~] = meshRead(fileName);
            aaa = load(outFileM);   % 'Vout', 'Nout', 'darFrames', 'partIDs';   % THESE ARE PARTS OF THE PREVIOUS LAYER
            Vout = aaa.Vout;
            Nout = aaa.Nout;
            darFrames = aaa.darFrames;
            partIDs = aaa.partIDs;
%             load(outFileAP);  % 'VAll', 'NAll', 'numPoints', 'areCentral' % THESE ARE ADDITIONAL POINTs generated from the mesh

            if size(Vout, 1) == 3
                Vout = Vout';
                Nout = Nout';
            end

            fring = compute_face_ring(F);
                
            nEl = size(Vout, 1);

            for j = 1:nEl 

                central_pos = Vout(j, :);
                NormalCentral = Nout(j, :);
                if layerID == 4
                    disp('TODO:define the central element!')
                end


                if layerID == 3 || layerID == 4        
                    %% compute neighbours
                    
                    % Create a list of faces
                    [facesId, nf] = ExtractAdjacentFaces(partIDs(j,2), nf, fring);
                    
                    % points normals and parts in the local neibourhood
                    points = zeros(1,3);
                    normals = zeros(1,3);
                    parts = zeros(1,1);
                    
                    % extract all parts that belong to these faces
                    startPoint = 1;
                    
                    for jj = 1:nf
                        ids = find(partIDs(:,2) == facesId(jj));
                        endPoint = startPoint + length(ids) - 1;
                        
                        points(startPoint:endPoint, :) = Vout(ids, :);
                        normals(startPoint:endPoint, :) = Nout(ids, :);
                        parts(startPoint:endPoint) = partIDs(ids, 1);
                        startPoint = endPoint + 1;
                    end

                    
                    if endPoint < 3  % not enought points in the neighbourhood
                        continue;
                    end
                    
                    % local frame of reference at this point
                    DF = darFrames(j, :);
                    DF = reshape(DF, [3,3]);

                    if layerID == 3
                        curDisp = DF(:,2);
                    elseif layerID == 4
                        curDisp = DF(:,3);
                    end
                end

                curOffset = zeros(3, 2);
                curOffset(:, 1) =   curDisp * offsetsConventional{layerID}; 
                curOffset(:, 2) = - curDisp * offsetsConventional{layerID};
                
                % find the closest point from  central_pos + curOffset(:, k)
                is_ok = true;
                
                for k = 1:2  % left and right neigbours

                    cent_cur = central_pos + curOffset(:, k)';
                    distances = ComputeEuclDists(cent_cur, points, 3);
                    % take the point with the smallest distance
                    minD = min(distances);
                    if minD > offsetsConventional{layerID}/3;
                        is_ok = false; 
                        break;
                    end
                    
                    inds = find(distances == minD);
                    if k == 1;
                        lefts = parts(inds(1));
                        V_left = points(inds(1), :);
                        NormalLeft = normals(inds(1),:);
                    else
                        rights = parts(inds(1));
                        V_right = points(inds(1), :);
                        NormalRight = normals(inds(1),:);
                    end
                    if minD > offsetsConventional{layerID}/2;
                        is_ok = false;
                    end                                      
                end
                
                if ~is_ok
                    continue;
                end

                if layerID == 3

                    if length(lefts) > 1  % find the one with the largest scale and shortest distance 
                        
                       disp('Something went wrong');
                       lenE = length(lefts);
%                        score = double(options.weight * lefts) - distsLeft;
                       id = find(score == max(score));
                       if length(id) > 1
                           id = id(1);
                       end
                       [lefts, indsXLeft, indsYLeft, distsLeft] = myFilter4(lefts,indsXLeft,indsYLeft,distsLeft, id);
                       depthsLeft= depthsLeft(id);
                    end

                    if length(rights) > 1  % find the one with the largest scale and shortest distance 
                       disp('Something went wrong');
                       lenE = length(rights);
%                        score = double(options.weight * rights) - distsRight;
                       id = find(score == max(score));
                       if length(id) > 1
                           id = id(1);
                       end
                       [rights, indsXRight, indsYRight, distsRight] = myFilter4(rights,indsXRight,indsYRight,distsRight, id);
                       depthsRight = depthsRight(id);
                    end

%                     if is_visualization
%                         XXX = [indsXLeft, indsXRight] - x + centre;
%                         YYY = [indsYLeft, indsYRight] - y + centre;
%                         ZZZ = [depthsLeft; depthsRight];
%                         scatter3(XXX, YYY, ZZZ, 'green');
%                         hold on
% 
%                     end

%                     if is_visualization
%                             plotFrame(NormalLeft, [indsXLeft-x+centre, indsYLeft-y+centre, depthsLeft], vecLen, vectColors, 1);
%                             plotFrame(NormalRight,[indsXRight-x+centre, indsYRight-y+centre, depthsRight], vecLen, vectColors, 1);
%                     end

                    % compute and quantize relative orienations of
                    % normals
                    Xtemp = DF(:,2);
                    Ytemp = DF(:,3);  % local x and y axis

                    % compute projections
                    angleXLeft =  90 - 180 * acos(dot(Xtemp,NormalLeft)/norm(NormalLeft)) /pi;
                    angleYLeft =  90 - 180 * acos(dot(Ytemp,NormalLeft)/norm(NormalLeft)) /pi;
                    angleXRight = 90 -  180 * acos(dot(Xtemp,NormalRight)/norm(NormalRight))/pi;
                    angleYRight = 90 - 180 * acos(dot(Ytemp,NormalRight)/norm(NormalRight))/pi;

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

                    %% compute relative depths of the pixels
                    T = eye(4); T(1,4) = -central_pos(1); T(2,4) = -central_pos(2); T(3,4) = -central_pos(3);
                    R = eye(4,4); R(1:3, 1) = Xtemp'; R(1:3, 2) = Ytemp'; R(1:3, 3) = NormalCentral';
%                     T2 = eye(4); T2(1,4) = central_pos(1); T2(2,4) = central_pos(2); T2(3,4) = central_pos(3);
                    M = inv(R)*T; 

                    % coordinates of the neighbours in the local frame
                    % of reference
                    coordsLeft = M * [V_left(1);  V_left(2);  V_left(3);  1];
                    coordsRight = M * [V_right(1);  V_right(2);  V_right(3); 1];

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

                    curCoords = [i, central_pos(1), central_pos(2), central_pos(3)];  % image, x, y
                    curCoords = uint16(curCoords);
                    outputCoords = [outputCoords; curCoords];
                    outputStatistics = [outputStatistics; line];
                    outputFrames = [outputFrames; DF(:)'];

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


%                     while length(facesId) < nf
%                         if length(facesId) == 1
%                             facesId = [facesId, fring{facesId}];
%                         else
%                             facesIdNew = facesId;
%                             for jj = 1:length(facesId)
%                                 facesIdNew = union(facesIdNew, fring{facesId(jj)});
%                             end
%                             facesId = facesIdNew;
%                         end
%                     end
%                     
%                     points = zeros(1,3);
%                     normals = zeros(1,3);
%                     
%                     % extract all parts that belong to these faces
%                     startPoint = 1;
%                     
%                     for jj = 1:nf
%                         ids = find(PI(:,2) == facesId(jj));
%                         endPoint = startPoint + length(ids) - 1;
%                         
%                         points(startPoint:endPoint, :) = V_el(ids, :);
%                         normals(startPoint:endPoint, :) = N(ids, :);
%                         startPoint = endPoint + 1;
%                     end




%                 % find a surface point that is closest to the [x,y,depthCentral]
%                 [xLeft, yLeft, zLeft] = findClosestPoint(I, [x + curOffset(1, 1), y + curOffset(2, 1),depthCentral + curOffset(3, 1)]);
%                 [xRight, yRight, zRight] = findClosestPoint(I, [x + curOffset(1, 2), y + curOffset(2, 2),depthCentral + curOffset(3, 2)]);

%                 if is_visualization
%                     shift = 30;
%                     centre = shift+1;
%                     II = I(y-30: y+30, x - 30: x+30);
%                     II(II < 1050) = NaN;
%                     % visualize this image in 3D
%                     surf(II, 'FaceColor',[0.3, 0.3 0.3], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceLighting', 'phong');
%                     camlight left
%                     axis equal;
%                     hold on
%                     curDepth = II(shift+1,shift+1);
% 
%                     vectColors = eye(3);
%                     vecLen = 20;
%                     plotFrame(V, [centre, centre, depthCentral], vecLen, vectColors, 3);
% 
%                     hold on
%                     plotFrame(eye(3), [centre, centre, depthCentral], vecLen, vectColors, 3);
% 
%                     hold on
%                     XXX = xs + centre;
%                     YYY = ys + centre;
%                     ZZZ = zs + depthCentral;
%                     scatter3(XXX, YYY, ZZZ, 'red');
%                     hold on
% 
%                 end
















