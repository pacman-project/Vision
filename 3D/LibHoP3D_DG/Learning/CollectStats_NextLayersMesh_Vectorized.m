% this is to learn co-occurrence statistics from meshes and point clouds


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [outputStatisticsAll, outputCoordsAll, outputFramesAll] = CollectStats_NextLayersMesh_PC(list_els, list_input, lenF, ...
                                                        nPrevClusters, layerID, depthStep, dataSetNumber, ...
                                                        maxRelDepth, cluster1Bounds, nClusters, offsetsConventional)
                                                 
    is_visualization = false;
    multX = 0.33;
    multY = 0.22;
    
    xyzStep = 0.001;
    angleStep = 0.01;
                                            
    disp('collecting co-occurrence statistics...');

    emptyCellID = nPrevClusters + 1; 
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end

    kkMax = 5;
    
    % numInitial = 10^4 * lenF;
    outputStatisticsAll = {};%zeros(numInitial, 27);
    outputCoordsAll = {};%zeros(numInitial, 4);
    outputFramesAll = {};%zeros(numInitial, 9);
    curTSAll = 0;
      

    for i = 1:lenF 
        
        curTS = 0;
        numEl = 0;

        fileName = list_els{i};
        outFileM = [fileName(1:end-4), '.mat'];   

        [~, F, ~] = meshRead(list_input{i});
        aaa = load(outFileM);               % 'Vout', 'Nout', 'likelihoods', 'darFrames', 'partIDs', 'NeiOut', 'pointIndexing', 'pointIDx' 
        Vout = aaa.Vout;
        Nout = aaa.Nout;
        darFrames = aaa.darFrames;
        partIDs = aaa.partIDs;
        pointIDx = aaa.pointIDx;
        
        
        clear('aaa');
%           load(outFileAP);  % 'VAll', 'NAll', 'numPoints', 'areCentral' % THESE ARE ADDITIONAL POINTs generated from the mesh


%             if is_visualization
%                 radii = 0.005 * ones(size(Vout, 2),1);
%                 [V_vis, F_vis] = prepareForVisualization(radii, Vout, Nout);
%                 VisualizeTriangulation(F_vis, V_vis);
%                 hold on
%             end

        if size(Vout, 1) == 3
            Vout = Vout';
            Nout = Nout';
        end

        nEl = size(Vout, 1);
        lenF = size(F, 2);
        fringAll = ComputeFringDeep(F, kkMax, lenF);   % to address: fringAll{2}{13}

        
%     outputStatisticsAll = zeros(numInitial, 27);
%     outputCoordsAll = zeros(numInitial, 4);
%     outputFramesAll = zeros(numInitial, 9);
        
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
        
        central_pos = zeros(nEl, 3);
        NormalCentral = zeros(nEl,3);
        centralInf = zeos(nEl, 2);
        M_page = zeros(4,4,nEl);
        
        points_page = 

        parfor j = 1:nEl   

            central_pos = Vout(j, :);
            NormalCentral = Nout(j, :);
            centralInf = partIDs(j, :);
            
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

            DF = darFrames(j, :);
            DF = reshape(DF, [3,3]);
            
            % compute and quantize relative orienations of normals
            Xtemp = DF(:,2);
            Ytemp = DF(:,3);  % local x and y axis

            %% compute relative depths of the pixels
            T = eye(4); T(1:3,4) = -central_pos';   % T(2,4) = -central_pos(2); T(3,4) = -central_pos(3);
            R = eye(4,4); R(1:3, 1) = Xtemp'; R(1:3, 2) = Ytemp'; R(1:3, 3) = NormalCentral';
%               T2 = eye(4); T2(1,4) = central_pos(1); T2(2,4) = central_pos(2); T2(3,4) = central_pos(3);

%             M = inv(R)*T;
            M = R\T;   
            M_page(:,:,j) = M;

        end









            % coordinates of the neighbours in the local frame
            % of reference
            coordsLeft =  M * [V_left, ones(size(V_left,1),  1)]';
            coordsRight = M * [V_right, ones(size(V_right,1),1)]';
            
            % filter out those that are far form the predicted positions
            idsTrueLeft  = abs(coordsLeft(1,:) - offsetsConventional{layerID}) < offsetsConventional{layerID} *multX  & abs(coordsLeft(2,:)) < offsetsConventional{layerID}*multY;
            idsTrueRight = abs(coordsRight(1,:) + offsetsConventional{layerID}) < offsetsConventional{layerID}*multX & abs(coordsRight(2,:)) < offsetsConventional{layerID}*multY;
            numLeft = nnz(idsTrueLeft);  numRight = nnz(idsTrueRight);
            
            if  numLeft ~= 0 && numRight ~= 0 
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

            
%     outputStatisticsAll = zeros(numInitial, 27);
%     outputCoordsAll = zeros(numInitial, 4);
%     outputFramesAll = zeros(numInitial, 9);
            lenVec = length(vec1);
            
            outputStatisticsTemp = zeros(27, lenVec);
            outputCoordsTemp = zeros(4, lenVec);
            outputFramesTemp = zeros(9, lenVec);
            
%     xyzStep = 0.001;
%     angleStep = 0.01;
            
            for k = 1:lenVec
                line = zeros(27, 1);
                
                line(1) = centralInf(1); 
                line(2) = leftInf(vec1(k), 1);     %leftInf = [part_number, fid, pointID]
                line(3:5) =   round(coordsLeft(1:3, vec1(k)) / xyzStep);
                NL = NormalLeft(vec1(k),:);
                line(6:8) = round([NL*Xtemp; NL*Ytemp; NL*NormalCentral']/angleStep);
                
                line(15) = rightInf(vec2(k), 1);     %leftInf = [part_number, fid, pointID]
                line(16:18) = round(coordsRight(1:3, vec2(k)) / xyzStep);
                NR = NormalRight(vec2(k),:);
                line(19:21) = round([NR*Xtemp; NR*Ytemp; NR*NormalCentral']/angleStep);

                curCoords = [i; pointIDx{centralInf(2)}(centralInf(3)); pointIDx{leftInf(vec1(k), 2)}(leftInf(vec1(k),3)); pointIDx{rightInf(vec2(k),2)}(rightInf(vec2(k),3))];  % meshID, poindIDx_c, poindIDx_l, poindIDx_r
                
                outputStatisticsTemp(:, k) = line;
                outputCoordsTemp(:, k) = curCoords;
                outputFramesTemp(:, k) = DF(:);
                
%                 outputCoords = [outputCoords, curCoords'];
%                 outputStatistics = [outputStatistics, line'];
%                 outputFrames = [outputFrames, DF(:)];
%                 curTS = curTS + 1;
            end
            outputStatistics{j} = int8(outputStatisticsTemp);
            outputCoords{j} = int32(outputCoordsTemp);
            outputFrames{j} = single(outputFramesTemp);

            if mod(j, 100) == 0
                str = [num2str(j), ' out of ', num2str(nEl)];
                disp(str);
            end
        end
        
        outputStatisticsAll{i} = [outputStatistics{:}]; % vertcat
        outputCoordsAll{i} = [outputCoords{:}];
        outputFramesAll{i} = [outputFrames{:}];
        
        if mod(i, 20) == 0
            str = ['Temp/stat',num2str(i), '.mat'];
            save(str, 'outputStatisticsAll', 'outputCoordsAll', 'outputFramesAll', '-v7.3');
            outputStatisticsAll = {};
            outputCoordsAll = {};
            outputFramesAll = {};
        end

%         outputStatisticsAll(curTSAll+1:curTSAll + curTS, :) = outputStatistics;
%         outputCoordsAll(curTSAll+1:curTSAll + curTS, :) = outputCoords;
%         outputFramesAll(curTSAll+1:curTSAll + curTS, :) = outputFrames;
%         curTSAll = curTSAll + curTS; 
        disp(['Image - ', num2str(i)]);
    end

%     outputStatisticsAll = outputStatisticsAll';
%     outputCoordsAll = outputCoordsAll';
%     outputFramesAll = outputFramesAll';
    
    disp('Statistics collection time is ...');
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



%             if is_visualization
%                 radii = 0.005 * ones(size(points, 1),1);
%                 [V_vis, F_vis] = prepareForVisualization(radii, points, normals);
%                 VisualizeTriangulation(F_vis, V_vis, [1.0, 0,0], [0.2, 0,0], 0.4);
%                 hold on
%             end

            % local frame of reference at this point

%             if is_visualization
%                 plotVector(curOffset(:, 1), 1.0, [1,0,0], central_pos);
%                 plotVector(curOffset(:, 2), 1.0, [0,1,0], central_pos);
%             end

%             if is_visualization
%                 XXX = [indsXLeft, indsXRight] - x + centre;
%                 YYY = [indsYLeft, indsYRight] - y + centre;
%                 ZZZ = [depthsLeft; depthsRight];
%                 scatter3(XXX, YYY, ZZZ, 'green');
%                 hold on
%                 plotFrame(NormalLeft, [indsXLeft-x+centre, indsYLeft-y+centre, depthsLeft], vecLen, vectColors, 1);
%                 plotFrame(NormalRight,[indsXRight-x+centre, indsYRight-y+centre, depthsRight], vecLen, vectColors, 1);
%             end



















