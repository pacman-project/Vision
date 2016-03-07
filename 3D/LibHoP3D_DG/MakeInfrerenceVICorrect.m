function  MakeInfrerenceVICorrect(list_els, list_input, lenFiles, ...
                                                        nPrevClusters, layerID, dataSetNumber, receptiveFieldRad, offsetsConventional, statMapProperties)
                                                    
    is_visualization = 0;
                                 
    str = ['Making inference of the layer ', num2str(layerID)];
    disp(str);
    
    strParts = ['Temp/','layer', num2str(layerID), '/partSelection', num2str(layerID), '.mat'];
    pp = load(strParts);      % 'partsOut', 'pairsAll', 'nNClusters', 'coverageOut'
    partsOut = pp.partsOut;
    pairsAll = pp.pairsAll;
    nNClusters = pp.nNClusters{layerID};
    
    numPairs = size(pairsAll, 1);
    numParts = size(partsOut, 1);
    
    % make triples structure
    triples = zeros(numPairs, numPairs, numPairs);
    for i = 1:numParts
        triples(partsOut(i, 1), partsOut(i, 2), partsOut(i, 3)) = i;
    end
    
%     aa = load('Temp/structCS.mat');
%     structCS = aa.structCS;
    
    strLayer = ['Temp/OR_node_layer_', num2str(layerID)];
    dd = load(strLayer);
    ORTable = dd.ORTable;
    
%   re-arrange list of pairs
    for i = 1:nPrevClusters
        for j = 1:nPrevClusters
            pairsAllCell{i,j} = [];
            pairIDsCell{i,j} = [];
        end
    end
    
    for i = 1:size(pairsAll,1)
        pairsAllCell{pairsAll(i, 1),pairsAll(i,2)} = [pairsAllCell{pairsAll(i, 1), pairsAll(i,2)}; pairsAll(i,:)];
        pairIDsCell{pairsAll(i, 1),pairsAll(i,2)} = [pairIDsCell{pairsAll(i, 1),pairsAll(i,2)}; i];
    end


    if layerID <= 4
        curvThresh = 0.5;
    elseif layerID <= 6
        curvThresh = -0.5;
    else
        curvThresh = -0.9;
    end
    
    multX = statMapProperties.multX * 1.8;
    multY = statMapProperties.multY * 1.8;
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
    end
  
    for i = 1:36 %37:lenFiles 
%
        fileName = list_els{i};
        if layerID == 3
            sc = 1;
            strScale = ['scale_', num2str(sc)];
            inFileM =  [fileName(1:end-4), strScale, '.mat'];
        else
            strScale = [];
            inFileM = [fileName(1:end-4), '.mat'];
        end
        inFilePS = [fileName(1:end-4), strScale, 'PS.mat'];
        
        strRep =   ['layer', num2str(layerID-1)];
        strRepTo = ['layer', num2str(layerID)];
        
        fileNameNext = strrep(fileName, strRep, strRepTo);
        outFileM =  [fileNameNext(1:end-4), '.mat'];    % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
%         outFilePS = [fileNameNext(1:end-4), 'PS.mat'];  % 'NeiOut'

        % load files
        
        if ~exist(inFileM, 'file')
            continue;
        end
        [V, F, ~] = meshRead(list_input{i});
        aaa = load(inFileM);
        Vout = aaa.Vout;
        Nout = aaa.Nout;
        darFrames = aaa.darFrames;
        partIDs = aaa.partIDs;
        partIDs(:, 5) = (1:size(partIDs, 1))';  % add realization ID in the end
        clear('aaa');
        
%         aaa = load(inFilePS);
%         NeiOut = aaa.NeiOut;
%         clear('aaa');
        
        VoutOut = {};
        NoutOut = {};
        darFramesOut = {};
        partIDsOut = {};
        triplesNei = {};
        
        % approximate position of each face centre
        V1s = V(:,F(1, :)); V2s = V(:,F(1, :)); V3s = V(:,F(1, :));
        Fpos = (V1s + V2s + V3s)/3;
        

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
        areFacesActive = zeros(1, size(F,2));

        for mm = 1:size(F,2)     % for each face shows which parts belong to it
            tempF2part{mm} = [];
        end
        
        for mm = 1:size(partIDs,1)
            tempF2part{partIDs(mm, 2)} = [tempF2part{partIDs(mm, 2)}, mm];
            areFacesActive(partIDs(mm, 2)) = 1;
        end
        areFacesActive = find(areFacesActive);
        lenFActive = length(areFacesActive);
        tempStruct = zeros(1, lenF);

        parfor iii = 1:lenFActive 
            
            VoutOutTemp = zeros(3, 100);
            NoutOutTemp = zeros(3, 100);
            darFramesOutTemp = zeros(100, 9);
            partIDsOutTemp = zeros(100, 5);
            triplesNeiTemp = zeros(3, 100);
            
            numOut = 0;
            curF = areFacesActive(iii);
                      
            % define distances in the GEODESIC RECEPTIVE FIELD   
            facesId = [curF, fringAll{kkMax}{curF}];

            % delete some of faces that are apparently too far away
            distsF = pdist2(Fpos(:, curF)', Fpos(:,facesId)');  % euclidian distances from this face distances
            facesId = facesId(distsF <= 2.2 * receptiveFieldRad);
            
            faces = F(:, facesId);
            faces1 = faces(:);
            [vertexAll,~,ic] = unique(faces1);
            Vtemp = V(:, vertexAll);
            faces2 = reshape(ic, size(faces));
            
            U = geodesicDistances(Vtemp, faces2, ic(1:3), nIterMax);
% % 
%             figure;
%             clf;
%             options.face_vertex_color = max(0, (receptiveFieldRad - U));
%             plot_mesh(Vtemp, faces2,options);
%             colormap jet(256);
%             lighting none;
%             axis equal;
% %             
            idsRF = U <= 1.2 * receptiveFieldRad;
            faces2T = sum(idsRF(faces2), 1);
            facesId = facesId(faces2T>0);
            
            idsCurFace = tempF2part{curF}';      % part ids for all parts at this face
            idsCurRF = [tempF2part{facesId}]';   % part ids for all parts at this receptive field
            
            points = Vout(idsCurRF, :);
            normals = Nout(idsCurRF, :);
            parts = partIDs(idsCurRF, :);   % part_ID, fid, pointID, pointIDx, realizationID
            frames = darFrames(idsCurRF, :);
            if length(idsCurRF) < 3  % not enought points in the neighbourhood
                continue;
            end
            
            for j = 1:length(idsCurFace) % for all parts in this face
                
                central_pos =   Vout(idsCurFace(j), :);
                NormalCentral = Nout(idsCurFace(j), :);
                centralInf = partIDs(idsCurFace(j), :);

                % local frame of reference at this point
                DFcur = darFrames(idsCurFace(j), :);
                DF = reshape(DFcur, [3,3]);
                
%                 figure
%                 VisualizeTriangulation(F(:, facesId), V)
%                 hold on
%                 plotFrame(DF, central_pos, receptiveFieldRad, eye(3), 3);
                
                cent_cur = central_pos;
                distances = sqrt(sum(([points(:,1)-cent_cur(1),points(:,2)-cent_cur(2),points(:,3)-cent_cur(3)]).^2, 2));  % Should USE geodesic distances here!!!
                inds = distances < 1.0 * receptiveFieldRad & distances > 0.5 * receptiveFieldRad;

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
%                 hold on 
%                 scatter3(V_all(:, 1), V_all(:, 2), V_all(:,3));
%                 a = 2;
%                 scatter3(points(:, 1), points(:, 2), points(:,3));
          
                
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
%                 scatter3(central_pos(1), central_pos(2), central_pos(3));
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

                % coordinates of the neighbours in the local frame of reference
                coordsAll =  M * [V_all, ones(size(V_all,1),  1)]';
                  
                % normals and frames in the local frame of reference
                NormalLocalAll = NormalAll * [Xtemp, Ytemp, NormalCentral']; %[NL*Xtemp; NL*Ytemp; NL*NormalCentral'] [NR*Xtemp; NR*Ytemp; NR*NormalCentral']
                
                framesLocalAll = zeros(size(framesAll));
                framesLocalAll(:, 1:3) = NormalLocalAll;
                framesLocalAll(:, 4:6) = framesAll(:, 4:6) * [Xtemp, Ytemp, NormalCentral'];
                framesLocalAll(:, 7:9) = framesAll(:, 7:9) * [Xtemp, Ytemp, NormalCentral'];

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
                NormalLeft = NormalLocalAll(idsLeft, :);
                NormalRight = NormalLocalAll(idsRight, :);
                framesLeft  = framesLocalAll(idsLeft, :);
                framesRight = framesLocalAll(idsRight, :);

                % filter out those that are far form the predicted positions
                if layerID == 3 
                    idsTrueLeft  = abs(coordsLeft(1,:) - offsetsConventional) < offsetsConventional *multX  & abs(coordsLeft(2,:)) < offsetsConventional*multY;
                    idsTrueRight = abs(coordsRight(1,:) + offsetsConventional) < offsetsConventional*multX & abs(coordsRight(2,:)) < offsetsConventional*multY;
                elseif layerID == 4 
                    idsTrueLeft  = abs(coordsLeft(1,:)) < offsetsConventional *multX  & abs(coordsLeft(2,:) - offsetsConventional) < offsetsConventional*multY;
                    idsTrueRight = abs(coordsRight(1,:)) < offsetsConventional*multX & abs(coordsRight(2,:) + offsetsConventional) < offsetsConventional*multY;
                elseif layerID == 5 || layerID == 7 % multX is interpreted as multD
                    idsTrueLeft  = abs(distancesLeft - offsetsConventional) <= offsetsConventional * multX & abs(coordsLeft(2,:))  < offsetsConventional*multY; 
                    idsTrueRight = abs(distancesRight- offsetsConventional) <= offsetsConventional * multX & abs(coordsRight(2,:)) < offsetsConventional*multY;
                elseif layerID == 6 || layerID == 8 % multY is interpreted as multD
                    idsTrueLeft  = abs(distancesLeft - offsetsConventional) <= offsetsConventional * multY & abs(coordsLeft(1,:))  < offsetsConventional*multX; 
                    idsTrueRight = abs(distancesRight- offsetsConventional) <= offsetsConventional * multY & abs(coordsRight(1,:)) < offsetsConventional*multX;
                end
%                 hold on
%                 scatter3(V_all(idsLeft(idsTrueLeft), 1), V_all(idsLeft(idsTrueLeft), 2), V_all(idsLeft(idsTrueLeft),3), 'filled', 'red');
%                 hold on
%                 scatter3(V_all(idsRight(idsTrueRight), 1), V_all(idsRight(idsTrueRight), 2), V_all(idsRight(idsTrueRight),3), 'filled', 'green');
%                 
%                 tt = framesAll(idsRight(idsTrueRight), :);
%                 numFrames = length(tt);
%                 centres = V_all(idsRight(idsTrueRight), :);
%                 for kkk = 1:numFrames
%                     V = reshape(tt(kkk, :), [3,3]);
%                     plotFrame(V, centres(kkk, :), receptiveFieldRad*0.3, eye(3), 3);
%                 end
                

                numLeft = nnz(idsTrueLeft);
                numRight = nnz(idsTrueRight);

                if  numLeft ~= 0 && numRight ~= 0
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
                
                %% matching!
                % compute score for each pair
                partIDsCur = [leftInf(:, 1); rightInf(:, 1)];
                
                if layerID <= 4  % convert normal into quaternions
                    
                    statistics = [NormalLeft; NormalRight];
                    normsV = sqrt(statistics(:,1).^2 + statistics(:,2).^2 + statistics(:,3).^2); % normalize the vectors
                    statistics = statistics./repmat(normsV, [1,3]);

                    % compute quaternion for each vector
                    statisticsOut = computeVectorQuaternionVectorized(statistics', [0,0,1]);
                    statisticsOut = qconj(statisticsOut, 1);
                    pointsCur = double([[coordsLeft(1:3, :)'; coordsRight(1:3, :)'], statisticsOut']);    %   all points
                    
                else %convert frame into quaternions
                    
                    statistics = [framesLeft; framesRight];
                    lenStat = size(statistics, 1);
                    normsV = sqrt(statistics(:,1).^2 + statistics(:,2).^2 + statistics(:,3).^2);
                    statistics(:,1:3) = statistics(:,1:3)./repmat(normsV, [1,3]);
                    normsX = sqrt(statistics(:,4).^2 + statistics(:,5).^2 + statistics(:,6).^2);
                    statistics(:,4:6) = statistics(:,4:6)./repmat(normsX, [1,3]);
                    normsY = sqrt(statistics(:,7).^2 + statistics(:,8).^2 + statistics(:,9).^2);
                    statistics(:,7:9) = statistics(:,7:9)./repmat(normsY, [1,3]);
                    
                    statisticsOut = zeros(lenStat, 4);
                    
                    for kk = 1:lenStat
                        dcm = [statistics(kk, 4:6); statistics(kk, 7:9); statistics(kk, 1:3)]; % [X; Y; N]  
                        q = dcm2q(dcm);
                        q = qnorm(q);
                        statisticsOut(kk, :) = q;
                    end
                    
                    pointsCur = double([[coordsLeft(1:3, :)'; coordsRight(1:3, :)'], statisticsOut]);
                end
                
                idsLRUnique = unique(partIDsCur);
                
                numPointsT = size(pointsCur, 1);
                D = 100 * ones(numPointsT, numPairs);
                
                for kk = 1:length(idsLRUnique)
                    idsTemp = partIDsCur == idsLRUnique(kk);
                    D(idsTemp, pairIDsCell{centralInf(1), idsLRUnique(kk)}) = Mixed_Eucl_Quat_cluster_dist(pointsCur(idsTemp, :), pairsAllCell{centralInf(1), idsLRUnique(kk)}(:, 3:end)); 
                end
                
                [leftM,  leftI] = min(D(1:numLeft, :),     [], 2);
                [rightM,rightI] = min(D(numLeft+1:end, :), [], 2);
                
                idsL = leftM <= 8.0;
                idsR = rightM <= 8.0;
                leftI = leftI(idsL); rightI = rightI(idsR);
                if isempty(leftI) || isempty(rightI)
                    continue;
                end
                NeiRealizIdsLeft = leftInf(idsL, 5);
                NeiRealizIdsRight = rightInf(idsR, 5);
                
                [leftM, ia_l] = unique(leftI); [rightM, ia_r] = unique(rightI);
                NeiRealizIdsLeft  = NeiRealizIdsLeft(ia_l);
                NeiRealizIdsRight = NeiRealizIdsRight(ia_r);
                
                vec1 = 1:length(leftM); vec2 = 1:length(rightM);
                [p,q] = meshgrid(vec1, vec2);
                vec1 = p(:);
                vec2 = q(:);
                partCands = zeros(1, 100);
                tripleCan = zeros(3, 100);
                numCands = 0;
                
                % do inference and matching
                for jjj = 1:length(vec2)
                    partCand = triples(leftM(vec1(jjj)), centralInf(1) , rightM(vec2(jjj)));
                    if partCand > 0
                        numCands = numCands + 1;
                        partCands(numCands) = partCand;
                        tripleCan(:, numCands) = [NeiRealizIdsLeft(vec1(jjj)), centralInf(5), NeiRealizIdsRight(vec2(jjj))];
                    end
                end
                
                if numCands >= 1
                    for tt = 1:min(numCands, 3) % save each part
                        numOut = numOut + 1;
                        VoutOutTemp(:, numOut) = central_pos';
                        NoutOutTemp(:, numOut) = NormalCentral';
                        triplesNeiTemp(:, numOut) = tripleCan(:, tt);
                        darFramesOutTemp(numOut, :) = DFcur;
                        partIDsOutTemp(numOut, :) = centralInf;
                        partIDsOutTemp(numOut, 1) = partCands(tt);
                    end
                end
            end
            if mod(iii, 200) == 0
                str = [num2str(iii), ' face out of ', num2str(lenFActive)];
                disp(str);
            end
            
            VoutOut{iii} = VoutOutTemp(:, 1:numOut);
            NoutOut{iii} = NoutOutTemp(:, 1:numOut);
            darFramesOut{iii} = darFramesOutTemp(1:numOut, :)';
            partIDsOut{iii} = partIDsOutTemp(1:numOut, :)';
            triplesNei{iii} = triplesNeiTemp(:, 1:numOut);
        end
        Vout = [VoutOut{:}];
        Nout = [NoutOut{:}];
        darFrames = [darFramesOut{:}]';
        partIDs = [partIDsOut{:}]';
        triplesNei = [triplesNei{:}];
        
        % compute NeiOutNext
        NeiOutNext = {};
        lenPartRealizationsNext = size(triplesNei, 2);
        
%         for ii = 1:lenPartRealizationsNext
%              NeiOutNext{ii} = [NeiOut{triplesNei(1, ii)}, NeiOut{triplesNei(2, ii)}, NeiOut{triplesNei(3, ii)}];
%         end
%         NeiOut = NeiOutNext;
        
        [folderName, ~] = getFolderName(outFileM);
        
        if ~exist(folderName, 'dir')
           mkdir(folderName);
        end
        
        
        % apply OR-Nodes
        partIDs(:, 1) = ORTable(partIDs(:, 1));
        partIDs = partIDs(:, 1:4);
        
        if layerID == 4 || layerID == 6
            partIDs(:, 1) = partIDs(:, 1) + 1;
        end
        
        

        
        
        
%         rr = load('Temp/faceIDsTrial.mat');
%         facesId = rr.facesId;
%         
%         color = structCS(:, facesId);
%         colorsLen = sqrt(color(1,:).^2 + color(2,:).^2 + color(3,:).^2);
%         color(1,:) = color(1,:)./colorsLen;
%         color(2,:) = color(2,:)./colorsLen;
%         color(3,:) = color(3,:)./colorsLen;
%         
%         trisurf(F(:, facesId)',V(1, :),V(2, :),V(3, :), 'FaceColor', 'flat', 'FaceVertexCData', color', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
%         hold on
%         
%         V_all = [];
%         Vout = Vout';
% 
%         for kk = 1:length(facesId)
%             ids = partIDs(:, 2) == facesId(kk); %& partIDs(:, 3) >= 2;
%             V_all = [V_all; Vout(ids, :)];
%         end
%         
%         scatter3(V_all(:, 1), V_all(:, 2), V_all(:,3)); 
%         
%         a = 2;
%         
% 
%         
%         Vout = Vout';
        
        
        
            
        if layerID == 4  
            sc = 2;
            strScale =  ['scale_', num2str(sc)];
            fileName =  list_els{i};
            inFileM1 =  [fileName(1:end-4), strScale, '.mat'];   % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
            inFilePS1 = [fileName(1:end-4), strScale, 'PS.mat']; % 'NeiOut'
            inFileNP1 =  [fileName(1:end-4), 'scale_1NP.mat'];

            inFileM1 = strrep(inFileM1, 'layer3', 'layer2');
            inFilePS1 = strrep(inFilePS1, 'layer3', 'layer2');
            inFileNP1 = strrep(inFileNP1, 'layer3', 'layer2');

%             tt = load(inFileNP1);
            aa = load(inFileM1);
            Vout = [Vout, aa.Vout];
            Nout = [Nout, aa.Nout];
            darFrames = [darFrames; aa.darFrames];
            partIDs = [partIDs; aa.partIDs];

%             aa = load(inFilePS1);
%             NO = aa.NeiOut;
%             pointScales = tt.pointScales;

%             % deleting points with scale == 1
%             for ii = 1:size(Vout, 2)
%                 if ii <= lenPartRealizationsNext
%                     temp = NeiOutNext{ii};
%                     scales = pointScales(temp);
%                     scales = scales > 1;
%                     temp = temp(scales);
%                     NeiOut{ii} = temp;
%                 else
%                     NeiOut{ii} = NO{ii - lenPartRealizationsNext};
%                 end
%             end  
        elseif layerID == 6
            sc = 3;
            strScale =  ['scale_', num2str(sc)];
            fileName =  list_els{i};
            inFileM1 =  [fileName(1:end-4), strScale, '.mat'];   % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
            inFilePS1 = [fileName(1:end-4), strScale, 'PS.mat']; % 'NeiOut'
            inFileNP1 =  [fileName(1:end-4), 'scale_1NP.mat'];

            inFileM1 = strrep(inFileM1,   'layer5', 'layer2');
            inFilePS1 = strrep(inFilePS1, 'layer5', 'layer2');
            inFileNP1 = strrep(inFileNP1, 'layer5', 'layer2');

%             tt = load(inFileNP1);
            aa = load(inFileM1);
            Vout = [Vout, aa.Vout];
            Nout = [Nout, aa.Nout];
            darFrames = [darFrames; aa.darFrames];
            partIDs = [partIDs(:, 1:4); aa.partIDs];

%             aa = load(inFilePS1);
%             NO = aa.NeiOut;
%             pointScales = tt.pointScales;
%             NeiOut = NeiOutNext;

% % %             for ii = 1:size(Vout, 2)
% % %                 if ii <= lenPartRealizationsNext
% % %                     temp = NeiOutNext{ii};
% % %                     scales = pointScales(temp);
% % %                     scales = scales > 2;
% % %                     temp = temp(scales);
% % %                     NeiOut{ii} = temp;
% % %                 else
% % %                     NeiOut{ii} = NO{ii - lenPartRealizationsNext};
% % %                 end
% % %             end  
        else
%             NeiOut = NeiOutNext;
        end
        
        
        
%         rr = load('Temp/faceIDsTrial.mat');
%         facesId = rr.facesId;
%         
%         color = structCS(:, facesId);
%         colorsLen = sqrt(color(1,:).^2 + color(2,:).^2 + color(3,:).^2);
%         color(1,:) = color(1,:)./colorsLen;
%         color(2,:) = color(2,:)./colorsLen;
%         color(3,:) = color(3,:)./colorsLen;
%         
%         trisurf(F(:, facesId)',V(1, :),V(2, :),V(3, :), 'FaceColor', 'flat', 'FaceVertexCData', color', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
%         hold on
%         
%         V_all = [];
%         Vout = Vout';
% 
%         for kk = 1:length(facesId)
%             ids = partIDs(:, 2) == facesId(kk); %& partIDs(:, 3) >= 2;
%             V_all = [V_all; Vout(ids, :)];
%         end
%         
%         scatter3(V_all(:, 1), V_all(:, 2), V_all(:,3)); 
%         
%         a = 2;
%         
% 
%         
%         Vout = Vout';
        
        
        
        scC{4} = 2;
        scC{6} = 3;
        if layerID == 4 || layerID == 6
            
            ids = partIDs(:, 3) >= scC{layerID};
            Vout = Vout(:, ids);
            Nout = Nout(:, ids);
            darFrames = darFrames(ids, :);
            partIDs = partIDs(ids, :);
        end
        
        
        if is_visualization
            
            vecLen = 0.1;
            vectColors = eye(3);
            
            numPartsVis = max(ORTable);
            if layerID == 4 || layerID == 6
                numPartsVis = numPartsVis + 1; % first do inference of those parts that are not 
            end
            for kk = 1:numPartsVis
                figure;
                iids = find(partIDs(:, 1) == kk);
                VV = Vout(:, iids);
                scatter3(VV(1,:), VV(2,:), VV(3,:)');

%                 hold on
%                 lenParts = length(iids);
%                 numNormals = min(100, lenParts);
%                 rr = randperm(lenParts, numNormals);
% 
%                 for kkkk = 1:numNormals
%                     VD = darFrames(rr(kkkk), :);
%                     VD = reshape(VD, [3,3]);
%                     Vcentre = Vout(:, rr(kkkk));
%                     plotFrame(VD, vecLen, vectColors, Vcentre);
%                     a = 2;
%                 end                
                a = 2;
            end
        end
        


        save(outFileM,  'Vout', 'Nout', 'darFrames', 'partIDs');    % needed for statistics collection
%         save(outFilePS, 'NeiOut');                                              % needed for part selection
        
    end
end


                
%                 ttt = unique(partIDsCur);
%                 for jj = 1:length(ttt)
%                     ids = partIDsCur == ttt(jj);
%                     pointsCurCur = pointsCur(ids, :);
%                 end



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










