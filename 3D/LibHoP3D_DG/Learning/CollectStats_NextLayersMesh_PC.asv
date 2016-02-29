% this is to learn co-occurrence statistics from meshes and point clouds


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

% [outputStatisticsAll, outputCoordsAll, outputFramesAll] =
function  CollectStats_NextLayersMesh_PC(list_els, list_input, lenFiles, ...
                                                        nPrevClusters, layerID, dataSetNumber, receptiveFieldRad, offsetsConventional, statMapProperties)
                                                 
%     is_visualization = false;

    numCombs = 25;
    disp('collecting co-occurrence statistics ... ');
    
    if layerID <= 4
        curvThresh = 0.5;
    elseif layerID <= 6
        curvThresh = -0.5;
    else
        curvThresh = -0.8;
    end
    
    multX = statMapProperties.multX;
    multY = statMapProperties.multY;
    xyzStep = statMapProperties.xyzStep;
    angleStep = statMapProperties.angleStep; 
                                                    

    emptyCellID = nPrevClusters + 1; 
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenFiles);
    end

    if layerID <= 4
        kkMax = 13;
        nIterMax = 100;
    elseif layerID <= 6 
        kkMax = 17;
        nIterMax = 150;
    elseif layerID <= 8 
        kkMax = 21;
        nIterMax = 220;
    elseif layerID <= 10 
        kkMax = 31;
        nIterMax = 300;
    end
    
%     outputStatisticsAll = {};
%     outputCoordsAll = {};
%     outputFramesAll = {};
    curTSAll = 0;
    
%     aa = load('Temp/structCS.mat');
%     structCS = aa.structCS;
      
    for i = 1:37%38:lenFiles 

        fileName = list_els{i};
        if layerID == 3
            sc = 1;
            strScale = ['scale_', num2str(sc)];
            outFileM =  [fileName(1:end-4), strScale, '.mat'];
        else
            outFileM = [fileName(1:end-4), '.mat'];
        end
        
        if ~exist(outFileM)
            continue;
        end

        % load files
        [V, F, ~] = meshRead(list_input{i});
        aaa = load(outFileM);
        Vout = aaa.Vout;
        Nout = aaa.Nout;
        darFrames = aaa.darFrames;
        partIDs = aaa.partIDs;
        
        
        partIDs(:, 5) = (1:size(partIDs, 1))';  % add realization ID in the end
        clear('aaa');

        if size(Vout, 1) == 3
            Vout = Vout';
            Nout = Nout';
        end
        if size(darFrames, 1) == 9
            darFrames = darFrames;
        end
        
        nEl = size(Vout, 1);
        lenF = size(F, 2);
        fringAll = ComputeFringDeep(F, kkMax, lenF);   % to address: fringAll{2}{13}

        outputStatistics = {};
        outputCoords = {};
        outputFrames = {};
        tempF2part = {};
        areFacesActive = zeros(1, size(F,2));

        for mm = 1:size(F,2) %for each face shows which parts belong to it
            tempF2part{mm} = [];
        end
        
        for mm = 1:size(partIDs,1)
            tempF2part{partIDs(mm, 2)} = [tempF2part{partIDs(mm, 2)}, mm];
            areFacesActive(partIDs(mm, 2)) = 1;
        end
        areFacesActive = find(areFacesActive);
        lenFActive = length(areFacesActive);
        
        % approximate position of each face centre
        V1s = V(:,F(1, :)); V2s = V(:,F(1, :)); V3s = V(:,F(1, :));
        Fpos = (V1s + V2s + V3s)/3;
        

        parfor iii = 1:lenFActive

            curF = areFacesActive(iii);
                      
            % define distances in the GEODESIC RECEPTIVE FIELD   
            facesId = [curF, fringAll{kkMax}{curF}];
            
            % delete some of faces that are apparently too far away
            distsF = pdist2(Fpos(:, curF)', Fpos(:,facesId)');  % euclidian distances from this face distances
            facesId = facesId(distsF <= 2.2 * receptiveFieldRad);

            faces = F(:, facesId);
            
%             figure
%             VisualizeTriangulation(faces, V);
%             hold on
            
            faces1 = faces(:);
            [vertexAll,~,ic] = unique(faces1);
            Vtemp = V(:, vertexAll);
            faces2 = reshape(ic, size(faces));
            
            U = geodesicDistances(Vtemp, faces2, ic(1:3), nIterMax);
            
%             figure;
%             clf;
%             options.face_vertex_color = max(0, (receptiveFieldRad - U));
%             plot_mesh(Vtemp, faces2,options);
%             colormap jet(256);
%             lighting none;
%             axis equal;
             
            idsRF = U <= 1.2 * receptiveFieldRad;
            faces2T = sum(idsRF(faces2), 1);
            facesId = facesId(faces2T>0);
%             
%             figure
%             VisualizeTriangulation(F(:, facesId), V);
%             hold on

            idsCurFace = tempF2part{curF}';      % part ids for all parts at this face
            idsCurRF = [tempF2part{facesId}]';   % part ids for all parts at this receptive field
            
            points = Vout(idsCurRF, :);
            normals = Nout(idsCurRF, :);
            parts = partIDs(idsCurRF, :);   % part_ID, fid, pointID, pointIDx, realizationID
            frames = darFrames(idsCurRF, :);
            if length(idsCurRF) < 3  % not enought points in the neighbourhood
                continue;
            end
            
            maxSizeStat = 2 * length(idsCurFace) * numCombs;
            outputStatisticsTemp = zeros(14, maxSizeStat);
            outputCoordsTemp = zeros(4, maxSizeStat);
            outputFramesTemp = zeros(10, maxSizeStat);
            curFrame = 0;
            curLine = 0;
            
            for j = 1:length(idsCurFace) % for all parts in this face
                
                central_pos =   Vout(idsCurFace(j), :);
                NormalCentral = Nout(idsCurFace(j), :);
                centralInf = partIDs(idsCurFace(j), :);

                % local frame of reference at this point
                DF = darFrames(idsCurFace(j), :);
                DF = reshape(DF, [3,3]);
                
%                 figure
% %                 VisualizeTriangulation(F(:, facesId), V);
%                 color = structCS(:, facesId);
%                 colorsLen = sqrt(color(1,:).^2 + color(2,:).^2 + color(3,:).^2);
%                 color(1,:) = color(1,:)./colorsLen;
%                 color(2,:) = color(2,:)./colorsLen;
%                 color(3,:) = color(3,:)./colorsLen;
% 
%                 trisurf(F(:, facesId)',V(1, :),V(2, :),V(3, :), 'FaceColor', 'flat', 'FaceVertexCData', color', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
%                 hold on
%                 plotFrame(DF, central_pos, receptiveFieldRad, eye(3), 3);
                
                                
                cent_cur = central_pos;
                distances = sqrt(sum(([points(:,1)-cent_cur(1),points(:,2)-cent_cur(2),points(:,3)-cent_cur(3)]).^2, 2));
                inds = distances < 1.0 * receptiveFieldRad & distances > 0.4 * receptiveFieldRad;

                if nnz(inds) < 3
                    continue;
                end

                % all points in the receptive field
                allInf = parts(inds, :);    % part_number, fid, pointID
                V_all = points(inds, :);  % vertex position
                
                framesAll = frames(inds, :);
                NormalAll = normals(inds,:);
                distances = distances(inds)';
%                 
%                 scatter3(points(:, 1), points(:, 2), points(:,3)); % all points from the adjacent faces
%                 hold on 
%                 scatter3(V_all(:, 1), V_all(:, 2), V_all(:,3)); % all points from the receptive field
%                 hold on
%                 a = 2;

% a = 2;
                
                
%                     if line(8)< -0.8 || line(21) < -0.8 
%                         disp(iii);
%                         
%                         hold on
%                         NormalCentral = NormalCentral*receptiveFieldRad;
%                         quiver3(central_pos(1), central_pos(2), central_pos(3), NormalCentral(1), NormalCentral(2), NormalCentral(3));
%                         
%                         hold on
%                         plotFrame(DF, central_pos, receptiveFieldRad, eye(3), 3);
%                         
%                         hold on
%                         quiver3(central_pos(1), central_pos(2), central_pos(3), NormalCentral(1), NormalCentral(2), NormalCentral(3));
%                         quiver3(central_pos(1), central_pos(2), central_pos(3), NormalCentral(1), NormalCentral(2), NormalCentral(3));
%                         
%                         a = 2;
%                     end
                
%                 scatter3(V_all(:,1), V_all(:,2), V_all(:,3), 'blue');
%                 hold on

%                 scatter3(central_pos(1), central_pos(2), central_pos(3), 'green');
%                 hold on
%                 plotFrame(DF, central_pos, receptiveFieldRad, eye(3), 3);
%                 axis equal

                % compute and quantize relative orienations of
                % normals
                if (NormalCentral * DF(:,1)) < 0.99
                    continue;
                end
                Xtemp = DF(:,2);
                Ytemp = DF(:,3);  % local x and y axis

                %% compute relative depths of the pixels
                T = eye(4); T(1:3,4) = -central_pos';   % T(2,4) = -central_pos(2); T(3,4) = -central_pos(3);
                R = eye(4,4); R(1:3, 1) = Xtemp'; R(1:3, 2) = Ytemp'; R(1:3, 3) = NormalCentral';

    %           T2 = eye(4); T2(1,4) = central_pos(1); T2(2,4) = central_pos(2); T2(3,4) = central_pos(3);
    %           M = inv(R)*T;
    
                M = R\T;

                % coordinates of the neighbours in the local frame
                % of reference
                coordsAll =  M * [V_all, ones(size(V_all,1),  1)]';

                zeroThresh = offsetsConventional * min(multY, multY);
                if layerID == 6
                    zeroThresh = zeroThresh/2;
                elseif layerID >= 7
                    zeroThresh = zeroThresh/6;
                end
                if mod(layerID, 2) == 1
                    idsLeft = coordsAll(1, :) >= zeroThresh;
                    idsRight = coordsAll(1, :) <= -zeroThresh;
                elseif mod(layerID, 2) == 0
                    idsLeft = coordsAll(2,:) >= zeroThresh;
                    idsRight = coordsAll(2, :) <= -zeroThresh;
                end

                coordsLeft = coordsAll(:, idsLeft);
                coordsRight = coordsAll(:, idsRight);
                distancesLeft  = distances(idsLeft);
                distancesRight = distances(idsRight);
                leftInf = allInf(idsLeft, :);
                rightInf = allInf(idsRight, :);
                NormalLeft = NormalAll(idsLeft, :);
                NormalRight = NormalAll(idsRight, :);
                framesLeft  = framesAll(idsLeft, :);
                framesRight = framesAll(idsRight, :);

                % filter out those that are far form the predicted positions
                if layerID == 3 
                    idsTrueLeft  = abs(coordsLeft(1,:) - offsetsConventional) < offsetsConventional *multX  & abs(coordsLeft(2,:)) < offsetsConventional*multY;
                    idsTrueRight = abs(coordsRight(1,:) + offsetsConventional) < offsetsConventional*multX & abs(coordsRight(2,:)) < offsetsConventional*multY;
                elseif layerID == 4 
                    idsTrueLeft  = abs(coordsLeft(1,:)) < offsetsConventional *multX  & abs(coordsLeft(2,:) - offsetsConventional) < offsetsConventional*multY;
                    idsTrueRight = abs(coordsRight(1,:)) < offsetsConventional*multX & abs(coordsRight(2,:) + offsetsConventional) < offsetsConventional*multY;
                elseif layerID == 5 || layerID == 7 || layerID == 9 % multX is interpreted as multD
                    idsTrueLeft  = abs(distancesLeft - offsetsConventional) <= offsetsConventional * multX & abs(coordsLeft(2,:))  < offsetsConventional*multY; 
                    idsTrueRight = abs(distancesRight- offsetsConventional) <= offsetsConventional * multX & abs(coordsRight(2,:)) < offsetsConventional*multY;
                elseif layerID == 6 || layerID == 8 || layerID == 10 % multY is interpreted as multD
                    idsTrueLeft  = abs(distancesLeft - offsetsConventional) <= offsetsConventional * multY & abs(coordsLeft(1,:))  < offsetsConventional*multX; 
                    idsTrueRight = abs(distancesRight- offsetsConventional) <= offsetsConventional * multY & abs(coordsRight(1,:)) < offsetsConventional*multX;
                end

                numLeft = nnz(idsTrueLeft);  
                numRight = nnz(idsTrueRight);

                if  numLeft ~= 0 && numRight ~= 0

                    if numLeft>numCombs
                        idsTrueLeft = find(idsTrueLeft);
                        idsCurRF = randperm(numLeft, numCombs);
                        idsTrueLeft = idsTrueLeft(idsCurRF);
                        numLeft = numCombs;
                    end
                    if numRight>numCombs
                        idsTrueRight = find(idsTrueRight);
                        idsCurRF = randperm(numRight, numCombs);
                        idsTrueRight = idsTrueRight(idsCurRF);
                        numRight = numCombs;
                    end
                    coordsLeft = coordsLeft(:, idsTrueLeft);
                    coordsRight = coordsRight(:, idsTrueRight);
                    leftInf = leftInf(idsTrueLeft', :);
                    rightInf = rightInf(idsTrueRight', :);
                    NormalLeft = NormalLeft(idsTrueLeft', :);
                    NormalRight = NormalRight(idsTrueRight', :);
                    framesLeft = framesLeft(idsTrueLeft', :);
                    framesRight = framesRight(idsTrueRight', :);
                else
                    continue;
                end
                
                % line shold have the following format:
                % partID_central 1, 
                % partID_left 2, 3-5: [x,y,z], 6:14 [N, X_axis, Y_axis] 
                % partID_right 15, 16-18: [x,y,z], 19-27: [N, X_axis, Y_axis]
                % X_axis, Y_axis might be unused for the first layers

                lenVec = numLeft + numRight;
                for k = 1:lenVec

                    curLine = curLine + 1;
                    line = zeros(14, 1); % zeros(27, 1);
                    line(1) = centralInf(1); 

                    if k<=numLeft
                        id = k;
                        line(2) = leftInf(id, 1);
                        line(3:5) = coordsLeft(1:3, id);  % round(coordsLeft(1:3, vec1(k)) / xyzStep);
                        NL = NormalLeft(id,:);
                        FL = framesLeft(id,:);
                        FLX = FL(4:6);
                        FLY = FL(7:9);
                        line(6:8)   = [NL*Xtemp; NL*Ytemp; NL*NormalCentral']; %round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                        line(9:11)  = [FLX*Xtemp; FLX*Ytemp; FLX*NormalCentral']; %round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                        line(12:14) = [FLY*Xtemp; FLY*Ytemp; FLY*NormalCentral']; %round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                        
                        curCoords = [i; centralInf(5); leftInf(id,5); 0];
                        if NL*NormalCentral' < curvThresh
                            curLine = curLine - 1;
                            continue;
                        end
                        if abs(line(10)) > 0.2 || abs(line(12)) > 0.2  % Gramma rule: to avoid gaps in the data
                            continue; 
                        end
                    else
                        id = k - numLeft;
                        line(2) = rightInf(id, 1);
                        line(3:5) = coordsRight(1:3, id); % round(coordsRight(1:3, vec2(k)) / xyzStep);
                        NR = NormalRight(id,:);
                        FR = framesRight(id,:);
                        FRX = FR(4:6);
                        FRY = FR(7:9);
                        line(6:8) = [NR*Xtemp; NR*Ytemp; NR*NormalCentral']; %round([NR*Xtemp; NR*Ytemp; NR*NormalCentral']/angleStep);
                        line(9:11)  = [FRX*Xtemp; FRX*Ytemp; FRX*NormalCentral']; %round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                        line(12:14) = [FRY*Xtemp; FRY*Ytemp; FRY*NormalCentral']; %round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                        
                        curCoords = [i; centralInf(5); 0; rightInf(id,5)];
                        if NR*NormalCentral' < curvThresh
                            curLine = curLine - 1;
                            continue;
                        end
                        if layerID == 6 || layerID == 8  % Gramma rule: to avoid gaps in the data
                            if abs(line(10)) > 0.2 || abs(line(12)) > 0.2
                                continue; 
                            end
                        end
                    end
%                     
%                     if line(5) > 0.004
%                         a = 2;
%                     end
%                      dcm = [line(9:11), line(12:14), line(6:8)]';
%                     line(3:5)

                    outputStatisticsTemp(:, curLine) = line;
                    outputCoordsTemp(:, curLine) = curCoords;
                end
                curFrame = curFrame + 1;
                outputFramesTemp(:, curFrame) = [centralInf(5); DF(:)];
            end
            
%             outputStatisticsTemp = outputStatisticsTemp(:, 1:curLine);
%             outputCoordsTemp = outputCoordsTemp(:, 1:curLine);
%             outputFramesTemp = outputFramesTemp(:, 1:curLine);
            
            outputStatistics{iii} = outputStatisticsTemp(:, 1:curLine);
            outputCoords{iii} = int32(outputCoordsTemp(:, 1:curLine));
            outputFrames{iii} = single(outputFramesTemp(:, 1:curFrame));
            
            if mod(iii, 500) == 0
                str = [num2str(iii), ' face out of ', num2str(lenFActive)];
                disp(str);
            end
        end
        
        outputCoordsAll = [outputCoords{:}];
        outputFramesAll = [outputFrames{:}];
        
        tempSt = [outputStatistics{:}]; % vertcat
        if ~isempty(tempSt)
            tempSt = convertWorldUnitsToFile(tempSt, xyzStep, angleStep);
        end
        outputStatisticsAll = tempSt;
        
        if mod(i, 1) == 0
            str = ['Temp/outputStatisticsAll_', num2str(i), '.mat'];
            save(str, 'outputStatisticsAll', '-v7.3');
            
            str = ['Temp/outputCoordsAll_', num2str(i) ,'.mat'];
            save(str, 'outputCoordsAll', '-v7.3');
            
%             outputStatisticsAll = {};
%             outputCoordsAll = {};
%             outputFramesAll = {};
        end

        disp(['Image - ', num2str(i)]);
    end
     
%     disp('Statistics collection time is ...');
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

function plotVector(V, vecLen, vectColors, point)

    XX = [point(1), point(1) + vecLen * V(1)];
    YY = [point(2), point(2) + vecLen * V(2)];
    ZZ = [point(3), point(3) + vecLen * V(3)];

    plot3(XX, YY, ZZ, 'Color', vectColors);
    hold on

end

function statistics = convertWorldUnitsToFile(statistics, xyzStep, angleStep)
    idsAngles = [(6:14)']; %(19:21)'
    idsCoords = [(3:5)']; %(16:18)'
    statistics(idsCoords, :) = round(statistics(idsCoords, :) / xyzStep);
    statistics(idsAngles, :) = round(statistics(idsAngles, :) / angleStep);
    statistics = int8(statistics);
end







%                 vec1 = 1:numLeft; vec2 = 1:numRight;
%                 [p,q] = meshgrid(vec1, vec2);
%                 vec1 = p(:);
%                 vec2 = q(:);
%                 lenVec = length(vec1);


%             curOffset = zeros(3, 2);
%             curOffset(:, 1) =   curDisp * offsetsConventional; 
%             curOffset(:, 2) = - curDisp * offsetsConventional;

%                 if is_visualization
%                     plotVector(curOffset(:, 1), 1.0, [1,0,0], central_pos);
%                     plotVector(curOffset(:, 2), 1.0, [0,1,0], central_pos);
%                 end

            % find the closest point from  central_pos + curOffset(:, k)



                % compute projections
%                 angleXLeft =  90 - 180 * acos(dot(Xtemp,NormalLeft)/norm(NormalLeft))   /pi;
%                 angleYLeft =  90 - 180 * acos(dot(Ytemp,NormalLeft)/norm(NormalLeft))   /pi;
%                 angleXRight = 90 - 180 * acos(dot(Xtemp,NormalRight)/norm(NormalRight)) /pi;
%                 angleYRight = 90 - 180 * acos(dot(Ytemp,NormalRight)/norm(NormalRight)) /pi;

%                 try
% 
%                 clusterXLeft = define1Cluster(angleXLeft, cluster1Bounds, nClusters);
%                 clusterYLeft  = define1Cluster(angleYLeft, cluster1Bounds, nClusters);
%                 clusterXRight = define1Cluster(angleXRight, cluster1Bounds, nClusters);
%                 clusterYRight = define1Cluster(angleYRight, cluster1Bounds, nClusters);
% 
%                 catch error1
%                     continue;
%                 end
% 
%                 if clusterXLeft <= 0 || clusterYLeft <= 0 || clusterXRight <=0 || clusterYRight <=0
%                     continue;
%                 end
% 
%                 % define left and right neighbours
%                 leftRO = compute2elementIndex(clusterXLeft, clusterYLeft, nClusters);    % relative orientation
%                 rightRO = compute2elementIndex(clusterXRight, clusterYRight, nClusters); % relative orientation



%                 offsetXLeft = abs(coordsLeft(1));
%                 offsetXRight = abs(coordsRight(1));
% 
%                 depthLeft = coordsLeft(3)  * offsetsConventional{layerID}/offsetXLeft;
%                 depthRight = coordsRight(3)* offsetsConventional{layerID}/offsetXRight;
% 
%                 if offsetXLeft == 0 || offsetXRight == 0;
%                     continue;
%                 end




%                     if length(lefts) > 1  % find the one with the shortest distance 
                        
%                        disp('Something went wrong');
%                        lenE = length(lefts);
% %                        score = double(options.weight * lefts) - distsLeft;
%                        id = find(score == max(score));
%                        if length(id) > 1
%                            id = id(1);
%                        end
%                        [lefts, indsXLeft, indsYLeft, distsLeft] = myFilter4(lefts,indsXLeft,indsYLeft,distsLeft, id);
%                        depthsLeft= depthsLeft(id);
%                     end
% 
%                     if length(rights) > 1  % find the one with the largest scale and shortest distance 
%                        disp('Something went wrong');
%                        lenE = length(rights);
% %                        score = double(options.weight * rights) - distsRight;
%                        id = find(score == max(score));
%                        if length(id) > 1
%                            id = id(1);
%                        end
%                        [rights, indsXRight, indsYRight, distsRight] = myFilter4(rights,indsXRight,indsYRight,distsRight, id);
%                        depthsRight = depthsRight(id);
%                     end






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



%                 radii = offsetsConventional/2 * ones(size(points, 1),1);
%                 [V_vis, F_vis, likelihoodsVis] = prepareForVisualization(radii, points, normals, 0.5*ones(size(radii)));
%                 VisualizeTriangulation(F_vis, V_vis, likelihoodsVis);
%                 hold on
%                 
%                 normals = normals / 100;
%                 for tt = 1:10:size(normals, 1)
%                     quiver3(points(tt, 1),points(tt, 2), points(tt, 3), normals(tt, 1), normals(tt, 2), normals(tt, 3));
%                     hold on
%                 end
%                 
%                 axis equal;
%                 a = 2;



                %distances1 = ComputeEuclDists(cent_cur, points, 3);  % TOO SLOW
%                distances1 = sqrt(sum((points - repmat(cent_cur, size(points,1), 1)).^2, 2));


%         points = Vout;
%         normals = Nout / 20;
%         for tt = 1:100:size(normals, 1)
%             quiver3(points(tt, 1),points(tt, 2), points(tt, 3), normals(tt, 1), normals(tt, 2), normals(tt, 3));
%             hold on
%         end



%             if is_visualization
%                 XXX = [indsXLeft, indsXRight] - x + centre;
%                 YYY = [indsYLeft, indsYRight] - y + centre;
%                 ZZZ = [depthsLeft; depthsRight];
%                 scatter3(XXX, YYY, ZZZ, 'green');
%                 hold on
%                 plotFrame(NormalLeft, [indsXLeft-x+centre, indsYLeft-y+centre, depthsLeft], vecLen, vectColors, 1);
%                 plotFrame(NormalRight,[indsXRight-x+centre, indsYRight-y+centre, depthsRight], vecLen, vectColors, 1);
%             end


%             figure
% %                 VisualizeTriangulation(F, V);
% %                 hold on
%             NormalCentral = NormalCentral/100;
%             quiver3(central_pos(1), central_pos(2), central_pos(3), NormalCentral(1), NormalCentral(2), NormalCentral(3));
%             hold on
%             VisualizeTriangulation(F(:, facesId), V);
% 
%             axis equal
% 
%             refine a list of faceIDs to measure topological distances!


