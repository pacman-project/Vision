% this is to learn co-occurrence statistics from meshes and point clouds


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatisticsAll, outputCoordsAll, outputFramesAll] = CollectStats_NextLayersMesh_PC_3D(list_els, list_input, lenF, ...
                                                        nPrevClusters, layerID, dataSetNumber, offsetsConventional, statMapProperties, curScale)
                                                 
%     is_visualization = false;

    numCombs = 6;
    disp('collecting co-occurrence statistics...');
    
%     stepI = 2;
    
    multX = statMapProperties.multX;
    multY = statMapProperties.multY;
    xyzStep = statMapProperties.xyzStep;
    angleStep = statMapProperties.angleStep; 
                                                    

    emptyCellID = nPrevClusters + 1; 
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end

    if layerID <= 4
        kkMax = 5;
    else 
        kkMax = 8;
    end
    outputStatisticsAll = {};
    outputCoordsAll = {};
    outputFramesAll = {};
    curTSAll = 0;
      
    for i = 1:lenF 

        fileName = list_els{i};
        
        if layerID == 3
            sc = 1;
            strScale = ['scale_', num2str(sc)];
            outFileM =  [fileName(1:end-4), strScale, '.mat'];
        end
        
%         outFileM = [fileName(1:end-4), '.mat'];   

        % load files
        [V, F, ~] = meshRead(list_input{i});
%         VisualizeTriangulation(F, V);
%         hold on
        
        aaa = load(outFileM);
        Vout = aaa.Vout;
        Nout = aaa.Nout;
        darFrames = aaa.darFrames;
        partIDs = aaa.partIDs;
        pointIDx = aaa.pointIDx;
        clear('aaa');

        if size(Vout, 1) == 3
            Vout = Vout';
            Nout = Nout';
        end
        
        nEl = size(Vout, 1);
        lenF = size(F, 2);
        fringAll = ComputeFringDeep(F, kkMax, lenF);   % to address: fringAll{2}{13}

        
        outputStatistics = {};
        outputCoords = {};
        outputFrames = {};
        tempF2part = {};

        for ii = 1:size(F,2)
            tempF2part{ii} = [];
        end
        for ii = 1:size(partIDs,1)
            tempF2part{partIDs(ii, 2)} = [tempF2part{partIDs(ii, 2)}, ii];
        end


        for j = 1:nEl   

            central_pos = Vout(j, :);
            NormalCentral = Nout(j, :);
            centralInf = partIDs(j, :);

            if layerID <= 6       
                %% compute neighbours

                facesId = fringAll{kkMax}{partIDs(j,2)};
                ids = [];
                if ~isempty(facesId)
%                   ids = ismember(partIDs(:,2), facesId); % too slow!
                    ids = [tempF2part{facesId}]';
                end
                
                if nnz(ids) < 3  % not enought points in the neighbourhood
                    continue;
                end  

                points = Vout(ids, :);
                normals = Nout(ids, :);
                parts = partIDs(ids, :);   % part_number, fid, pointID
                
                % local frame of reference at this point
                DF = darFrames(j, :);
                DF = reshape(DF, [3,3]);

                if layerID == 3 || layerID == 5
                    curDisp = DF(:,2);
                elseif layerID == 4 || layerID == 6
                    curDisp = DF(:,3);
                end
            end

            cent_cur = central_pos;
            distances = sqrt(sum(([points(:,1)-cent_cur(1),points(:,2)-cent_cur(2),points(:,3)-cent_cur(3)]).^2, 2));
    
            % take points that are not far away
            idsD = distances > 0.7 * offsetsConventional & distances < 0.3 * offsetsConventional;
            points = points(idsD, :);

            % compute and quantize relative orienations of
            % normals
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
            coordsLeft =  M * [V_left, ones(size(V_left,1),  1)]';
            coordsRight = M * [V_right, ones(size(V_right,1),1)]';
            
            % filter out those that are far form the predicted positions
            if layerID == 3 || layerID == 5
                idsTrueLeft  = abs(coordsLeft(1,:) - offsetsConventional) < offsetsConventional *multX  & abs(coordsLeft(2,:)) < offsetsConventional*multY;
                idsTrueRight = abs(coordsRight(1,:) + offsetsConventional) < offsetsConventional*multX & abs(coordsRight(2,:)) < offsetsConventional*multY;
            elseif layerID == 4 || layerID == 6
                idsTrueLeft  = abs(coordsLeft(1,:)) < offsetsConventional *multX  & abs(coordsLeft(2,:) - offsetsConventional) < offsetsConventional*multY;
                idsTrueRight = abs(coordsRight(1,:)) < offsetsConventional*multX & abs(coordsRight(2,:) + offsetsConventional) < offsetsConventional*multY;
            end
            numLeft = nnz(idsTrueLeft);  numRight = nnz(idsTrueRight);
            
            if  numLeft ~= 0 && numRight ~= 0
                
                if numLeft>numCombs
                    idsTrueLeft = find(idsTrueLeft);
                    ids = randperm(numLeft, numCombs);
                    idsTrueLeft = idsTrueLeft(ids);
                    numLeft = numCombs;
                end
                if numRight>numCombs
                    idsTrueRight = find(idsTrueRight);
                    ids = randperm(numRight, numCombs);
                    idsTrueRight = idsTrueRight(ids);
                    numRight = numCombs;
                end
                coordsLeft = coordsLeft(:, idsTrueLeft);
                coordsRight = coordsRight(:, idsTrueRight);
                leftInf = leftInf(idsTrueLeft', :);
                rightInf = rightInf(idsTrueRight', :);
                NormalLeft = NormalLeft(idsTrueLeft', :);
                NormalRight = NormalRight(idsTrueRight', :);
            else
                continue;
            end
            
            % make a list of pairs of lefts and rights
            vec1 = 1:numLeft; vec2 = 1:numRight;
            [p,q] = meshgrid(vec1, vec2);
            vec1 = p(:);
            vec2 = q(:);

            % line shold have the following format:
            % partID_central 1, 
            % partID_left 2, 3-5: [x,y,z], 6:14 [N, X_axis, Y_axis] 
            % partID_right 15, 16-18: [x,y,z], 19-27: [N, X_axis, Y_axis]
            % X_axis, Y_axis might be unused for the first layers

            
            lenVec = length(vec1);
            outputStatisticsTemp = zeros(27, lenVec);
            outputCoordsTemp = zeros(4, lenVec);
            outputFramesTemp = zeros(9, lenVec);
            
            for k = 1:lenVec
                line = zeros(27, 1);
                
                line(1) = centralInf(1); 
                line(2) = leftInf(vec1(k), 1);
                line(3:5) = coordsLeft(1:3, vec1(k));  % round(coordsLeft(1:3, vec1(k)) / xyzStep);
                NL = NormalLeft(vec1(k),:);
                line(6:8) = [NL*Xtemp; NL*Ytemp; NL*NormalCentral']; %round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                
                line(15) = rightInf(vec2(k), 1);
                line(16:18) = coordsRight(1:3, vec2(k)); % round(coordsRight(1:3, vec2(k)) / xyzStep);
                NR = NormalRight(vec2(k),:);
                line(19:21) = [NR*Xtemp; NR*Ytemp; NR*NormalCentral']; %round([NR*Xtemp; NR*Ytemp; NR*NormalCentral']/angleStep);

                curCoords = [i; pointIDx{centralInf(2)}(centralInf(3)); pointIDx{leftInf(vec1(k), 2)}(leftInf(vec1(k),3)); pointIDx{rightInf(vec2(k),2)}(rightInf(vec2(k),3))];  % meshID, poindIDx_c, poindIDx_l, poindIDx_r
                outputStatisticsTemp(:, k) = line;
                outputCoordsTemp(:, k) = curCoords;
                outputFramesTemp(:, k) = DF(:);
            end

            if mod(j, 10^4) == 0
                str = [num2str(j), ' out of ', num2str(nEl)];
                disp(str);
            end
            
            outputStatistics{j} = outputStatisticsTemp;
            outputCoords{j} = int32(outputCoordsTemp);
            outputFrames{j} = single(outputFramesTemp);
        end
        
        outputCoordsAll{i} = [outputCoords{:}];
        outputFramesAll{i} = [outputFrames{:}];
        
        tempSt = [outputStatistics{:}]; % vertcat
        tempSt = convertWorldUnitsToFile(tempSt, xyzStep, angleStep);
        outputStatisticsAll{i} = tempSt;
        
        if mod(i, 1) == 0
            str = ['Temp/outputStatisticsAll_', num2str(i), '.mat'];
            save(str, 'outputStatisticsAll', '-v7.3');
            
            str = ['Temp/outputCoordsAll_', num2str(i) ,'.mat'];
            save(str, 'outputCoordsAll', '-v7.3');
            
            outputStatisticsAll = {};
            outputCoordsAll = {};
            outputFramesAll = {};
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
    idsAngles = [(6:8)'; (19:21)'];
    idsCoords = [(3:5)'; (16:18)'];
    statistics(idsCoords, :) = round(statistics(idsCoords, :) / xyzStep);
    statistics(idsAngles, :) = round(statistics(idsAngles, :) / angleStep);
    statistics = int8(statistics);
end




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






%         points = Vout;
%         normals = Nout / 20;
%         for tt = 1:100:size(normals, 1)
%             quiver3(points(tt, 1),points(tt, 2), points(tt, 3), normals(tt, 1), normals(tt, 2), normals(tt, 3));
%             hold on
%         end








