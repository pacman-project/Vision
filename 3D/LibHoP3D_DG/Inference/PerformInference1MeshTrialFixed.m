% this is to perform inference of the first layer parts with FIXED size of
% the part

function PerformInference1MeshTrialFixed(list_input, lenFiles, patchRad, receptiveFieldRad, is_overwrite, strE) %, crossScaleStructure, curScale)
  
    is_visualization = true;
    
    kkMax = 16;
    maxNumPointDF = 600;

    scales = [1.0, 3.0, 9.0];
    frameOfReferenceMultiplier = [1.3, 1.0, 0.8];
    inferencePoints = [1,2,3];  % start inference at these points only
    
    a{1} = [1,32,38,40,43,51,57,72];  % bottles
    a{2} = [8 16,30,31];              % bowls
    a{3} = [33, 34, 35, 36, 37];      % box
    a{4} = [41, 42];                  % cup
    a{5} = [45, 46, 47, 48, 49, 50, 52];  % cutting board
    a{6} = [39];                      % can
    a{7} = [13, 14 15 44];            % tray
    a{8} = [2 3 4 73 74 75 76];       % plate 
    a{9} = [6,7,9];  % tea cup
    a{10} = [5];  % salt
    a{11} = [10, 11 12, 65];  % teapot
    a{12} = [17:29 60];  % vase
    a{13} = [53:56 58 59]; % frying pan
    a{14} = [61:64, 66]; % Jug
    a{15} = [67:71]; % mug
    
%     filesExist = zeros(3,75);
     
    for kk =  1:lenFiles
        
        
        fullFileName = list_input{kk};
        [V, F, ~] = meshRead(fullFileName);
%         figure
%         VisualizeTriangulation(F, V);


%         for sc = 3:-1:1
%             
%             [~, fileName] = getFolderName(fullFileName);
%             strScale = ['scale_', num2str(sc)];
%             outFileM =  [strE, '/', fileName(1:end-4), strScale, '.mat'];
%             
%             if exist(outFileM)
%                filesExist(sc, kk) = 1;
%             end 
%         end

%         if kk == 76
%             a =2;
%         end
%         continue;
        
        
        str = ['Inference from model ', num2str(kk)];
        disp(str);

        fullFileName = list_input{kk};        
        [V, F, ~] = meshRead(fullFileName);

        lenV = size(V, 2);
        lenF = size(F, 2);
        fringAll = ComputeFringDeep(F, kkMax, lenF);
        fring = fringAll{kkMax};
        
        fileNameNormals = [fullFileName(1:end-4), '_Normal.mat'];
        
        [N, ~] = compute_normal(V,F, fileNameNormals);
        
%         VisualizeTriangulation(F, V);
%         hold on  
%         for jj = 1:200:lenF
%             fcur = jj;
% %             quiver3(FposAll(1, jj), FposAll(2, jj), FposAll(3, jj), NF(1, jj), NF(2,jj), NF(3,jj));
%             
%             for jjj = 1:3
%                 vc = V(:, F(jjj, fcur));
%                 nc = N(:, F(jjj, fcur));                       
%                 quiver3(vc(1), vc(2), vc(3), nc(1), nc(2), nc(3));
%                 hold on
%             end
%             a = 2;
%         end
        
%         N = CheckNormals(V, F, N, NF);
%         N = N';
        
        %% make sure that size of the object is 1 at its larger dimension

        threshPlanar = [18, 13, 8]; % degrees
        gridStep = patchRad/3.0;

%         % visualize the mesh and some of the normals
%         [V,F] = check_face_vertex(V,F); 

        F1= F(1,:); F2 = F(2,:); F3 = F(3,:);
        V1 = V(:,F1); V2 = V(:,F2); V3 = V(:,F3);
        FposAll = (V1+V2+V3)/3;  % centre of each face
        
        crossScaleStructure = zeros(lenF, 3);
        crossScaleStructureDilated = zeros(lenF, 3);
        [VAll, NAll, pointIDx, numPoints, areCentral, areCentralAdd] = compute_AdditionalPointsRegular(V, F, N, gridStep);
        
        pointScales = [areCentral{:}]';
%         VisualizeTriangulation(F, V);

        for sc = 3:-1:1
            
            Vout = [];
            Nout = [];
            Fout = [];   % shows to which face does this part belong
            darFrames = [];
            likelihoods = [];
            partScales = [];
            NeiOutTemp = {};
            
            [~, fileName] = getFolderName(fullFileName);
            strScale = ['scale_', num2str(sc)];
            outFileM =  [strE, '/', fileName(1:end-4), strScale, '.mat'];
            outFilePS = [strE, '/', fileName(1:end-4), strScale, 'PS.mat'];
            outFileLH = [strE, '/', fileName(1:end-4), strScale, 'LH.mat'];
            outFileNP = [strE, '/', fileName(1:end-4), strScale, 'NP.mat'];
            
            patchRadTemp = patchRad * scales(sc);
            receptiveFieldRadTemp = receptiveFieldRad * scales(sc);
            
            vecLen = patchRadTemp * 4;
            vectColors = eye(3);
            curvThresh = 1/(20*patchRadTemp);
            
            str = ['Inference at scale ', num2str(sc)];
            disp(str);
            
            % dilate
            if sc == 2 
                crossScaleStructureDilated(:, sc+1) = dilateCrossScaleStructureTopologic(crossScaleStructure(:, sc+1), F, V, receptiveFieldRadTemp);
            elseif sc == 1
                temp = crossScaleStructure(:, sc+1) | crossScaleStructure(:, sc+2);
                crossScaleStructureDilated(:, sc+1) = dilateCrossScaleStructureTopologic(temp, F, V, 2 * receptiveFieldRadTemp);
            end
            
%             if sc == 1 || sc == 3
%                 areCentralCur = areCentral;
%             elseif sc == 2
%                 areCentralCur = areCentralAdd;
%             end

            areCentralCur = areCentral;
            
%             rr = load('Temp/faceIDsTrial.mat');
%             facesIdMy = rr.facesId;
                        
             parfor i = 1:lenF
                
%                 i = facesIdMy(jj);
                
                % check if something is inferred on the previous scales
                if sc <= 2
                    if crossScaleStructureDilated(i, sc+1) == 1
                        continue;
                    end
                end

                facesId = [i, fring{i}];
                Fpos = FposAll(:,facesId)';
                Fdist = pdist2(Fpos(1,:), Fpos);

                idsF = Fdist <= 3 * receptiveFieldRadTemp;
                facesId = facesId(idsF);

                points = zeros(3, 100);
                Normals = zeros(3, 100);
%                 FaceIDs = zeros(1, 100);
%                 pointIDs = zeros(1, 100);
                pointIDXs = zeros(1, 100);

                % EXTRACT POINTS AND NORMALS THAT BELONG TO THESE FACES    
                startPoint = 1;
                endPoint = 0;
                for j = 1:length(facesId)
                    if numPoints(facesId(j)) > 0
                        areCentralF = areCentralCur{facesId(j)};
                        idsCent = find(areCentralF >= inferencePoints(sc)); 
                        endPoint = startPoint + length(idsCent) - 1;
                        points(:,  startPoint:endPoint) = VAll{facesId(j)}(idsCent, :)';
                        Normals(:, startPoint:endPoint) = NAll{facesId(j)}(idsCent, :)';
%                         FaceIDs(startPoint:endPoint) = facesId(j);
%                         pointIDs(startPoint:endPoint) = idsCent;
                        pointIDXs(startPoint:endPoint) = pointIDx{facesId(j)}(idsCent);
                        startPoint = endPoint + 1;
                    end
                end
                if endPoint == 0
                    continue;
                end

                points = points(:, 1:endPoint);
                Normals = Normals(:, 1:endPoint);
%                 FaceIDs = FaceIDs(1:endPoint);
%                 pointIDs = pointIDs(1:endPoint);
                pointIDXs = pointIDXs(1:endPoint);
                
                areCentralF = areCentralCur{i};
                idsCent = find(areCentralF >= inferencePoints(sc));
                numM = nnz(idsCent);

                for k = 1:length(idsCent)

                    pointID = idsCent(k);
                    pointScale = areCentralF(pointID);
                    Vcentre = VAll{i}(idsCent(k), :)';
                    Ncentre = NAll{i}(idsCent(k), :)';

                    % COMPUTE DISTANCES FROM THE CENTRAL POINT      
                    dists = sqrt(sum((points - repmat(Vcentre, 1, size(points,2))).^2, 1));
                    % COMPUTE DISTANCES FROM THE CENTRAL NORMAL (via dot product)
                    dotProds = sum(Normals .* repmat(Ncentre, 1, size(points,2)), 1);

                    ids = find(dists>0 & dists <= patchRadTemp & dotProds > -0.2);
                    idsDF  = find(dists>0 & dists <= (frameOfReferenceMultiplier(sc) * receptiveFieldRadTemp) & dotProds > -0.2);

                    if length(ids) < 4 || length(idsDF) < 6
                        continue;
                    end

                    % this is done in order to speed up  

                    pointsTemp = points(:, ids);
                    NormalsTemp = Normals(:, ids);   

                    %% TRY TO SPEED UP

                    nc = repmat(Ncentre, [1, size(pointsTemp, 2)]);
                    an = abs( acos(sum(nc.*NormalsTemp,1))  * 180 / pi); % dot(nc, NormalsTemp, 1)
                    err = sum(an) / size(pointsTemp, 2);


                    if err < threshPlanar(sc)

                        likelihood = max(0, 1 - err/threshPlanar(sc));
                        Pind = pointIDx{i}(pointID); % global ID of this point
                        
                        Vout = [Vout, Vcentre];
                        Nout = [Nout, Ncentre];
                        Fout = [Fout, [i; pointID; Pind]];
                        partScales = [partScales; pointScale];
                        likelihoods = [likelihoods; likelihood];

%                         if length(idsDF) > maxNumPointDF
%                             idsTemp = randperm(length(idsDF), maxNumPointDF);
%                             idsDF = idsDF(idsTemp);
%                         end
                        
                        pointsDFTemp = points(:, idsDF);
                        normalsDFTemp = Normals(:, idsDF);
                        
%                         colormap jet;
%                         cmap = colormap;
%                         numColors = size(cmap, 1);
%                         distsDF = dists(idsDF);
%                         distMax = max(distsDF);
%                         distsDF = round(numColors * distsDF/distMax);
%                         distsDF(distsDF == 0) = 1;
%                         distsDF(distsDF > numColors) = numColors;
%                         colours = cmap(distsDF, :);
%                         sizeS = 2 * ones(length(distsDF),1);
%                         pointsDFTemp = points(:, idsDF);
%                         scatter3(pointsDFTemp(1,:), pointsDFTemp(2, :), pointsDFTemp(3,:), sizeS, colours);
%                         axis equal;
                        
                        pointsDFTemp = pointsDFTemp - repmat(Vcentre, [1, length(idsDF)]);
                        [VV, ~] = computeDarbouxFrame_V_N(Ncentre, pointsDFTemp, normalsDFTemp, curvThresh);
                        
                        sDF = sum(VV(:));
                        if isnan(sDF)
                            str = ['ERROR1:', num2str(i), ' ', num2str(k)];
                            disp(str);
                        end
                        if (VV(:, 1)' * Ncentre) < 0.99
                            str = ['ERROR2:', num2str(i), ' ', num2str(k)];
                            disp(str);
                        end

                        darFrames =  [darFrames,  VV(:)];

    %                     plotFrame(VV, vecLen, vectColors, Vcentre);
    %                     hold on

                        % get a list of points that belong to this planar patch
%                         FaceIDsTemp = FaceIDs(ids);
%                         pointIDsTemp = pointIDs(ids);
%                         for nn = 1:length(FaceIDsTemp)
%                             Pind = [Pind, pointIDx{FaceIDsTemp(nn)}(pointIDsTemp(nn))];
%                         end
                        pointIDXsTemp = pointIDXs(ids);
                        NeiOutTemp{i}{pointID} = uint32([Pind, pointIDXsTemp]);
                        %NeiOut{pointIDxCur} = Pind;
                    end

                end

                if mod(i, 500) == 0
                    str = [num2str(i), ' out of ', num2str(lenF)];
                    disp(str);
                end 
            end
            

            NeiOut = {};
            tic
            for i = 1:size(Fout, 2)
                temp = Fout(:, i);
                NeiOut{i} = NeiOutTemp{temp(1)}{temp(2)};
            end
            toc

            
%             lenNO = length(NeiOutTemp);
%             tic
%             for i = 1:size(pointIndexing, 1)
%                 temp = pointIndexing(i, :);
%                 if temp(1) < lenNO && temp(2) < length(NeiOutTemp{temp(1)})
%                     NeiOut{i} = NeiOutTemp{temp(1)}{temp(2)};
%                 else
%                     NeiOut{i} = [];
%                 end
%             end
%             toc
            
            darFrames = darFrames';
            Fout = Fout';

            % VISUALIZE THE RECONSTRUCTION 
            lenParts = length(likelihoods);
            if lenParts == 0
                continue;
            end
            
            radii = patchRadTemp * ones(size(likelihoods));
            
            crossScaleStructure(Fout(:,1), sc) = 1;

            if is_visualization
                figure
                cmap = load('settings/colormapGray.mat');  cmap = cmap.cmap;
                is_subset = false;
                
                if is_subset
                    iiddss = randperm(lenParts,50);
                    Vout = Vout(:, iiddss);
                    Nout = Nout(:, iiddss);
                    radii = radii(iiddss);
                    likelihoods = likelihoods(iiddss);
                end
                [V_vis, F_vis, likelihoodsOut] = prepareForVisualization(radii, Vout, Nout, likelihoods);
                VisualizeTriangulation(F_vis, V_vis, likelihoodsOut, cmap);
                hold on

                numNormals = min(200, lenParts);
                rr = randperm(lenParts, numNormals);

                for kkkk = 1:numNormals
                    VD = darFrames(rr(kkkk), :);
                    VD = reshape(VD, [3,3]);
                    Vcentre = Vout(:, rr(kkkk));
                    plotFrame(VD, vecLen, vectColors, Vcentre);
                    a = 2;
                end
                hold on

            end  % if is_visualization

            partIDs = ones(lenParts, 4);

            if ~isempty(partIDs)
                partIDs(:,2:4) = Fout;
            end

            if ~exist(strE, 'dir')
               mkdir(strE);
            end

            if is_overwrite || ~exist(outFileM, 'file')
                doWrite = 1;
            end
            tic
            if doWrite
                
                partIDs(:, 3) = partScales;
                Vout = single(Vout);
                Nout = single(Nout);
                darFrames = single(darFrames);
                partIDs = int32(partIDs);
                likelihoods = single(likelihoods);
                
                a = 2;
                
                save(outFileM,  'Vout', 'Nout', 'darFrames', 'partIDs', 'partScales');  % needed for statistics collection
%                 save(outFilePS, 'NeiOut');                     % needed for part selection
%                 save(outFileLH, 'likelihoods', 'partIDs');     % needed for  reviseFaces
%                 save(outFileNP, 'numPoints', 'pointScales');   % needed for  reviseFaces
            end
            toc
        end % end scale loop
        
        if mod(kk, 5) == 0
            a = 2;
        end
    end
end


function plotFrame(V, vecLen, vectColors, point)
    for ii = 1:3
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [point(1), point(1) + vecLen * curVect(1)];
        YY = [point(2), point(2) + vecLen * curVect(2)];
        ZZ = [point(3), point(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on
    end
end

function numPoints = computeNumPoints(areas, q30, q90, maxPoints)
    
    for i = 1:length(areas)
        if areas(i) < q30
            numPoints(i) = 0;
        elseif areas(i) > q90
            numPoints(i) = maxPoints;
        else
            numPoints(i) = round(maxPoints*(areas(i) - q30)/(q90 - q30));
        end
    end
end







%         aa = load('Temp/structCS.mat');
%         structCS = aa.structCS;    
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
%         Fout = Fout';
% 
%         for kk = 1:length(facesId)
%             ids = Fout(:, 1) == facesId(kk); %& partIDs(:, 3) >= 2;
%             V_all = [V_all; Vout(ids, :)];
%         end
%         
%         scatter3(V_all(:, 1), V_all(:, 2), V_all(:,3)); 
%         
%         a = 2;
%         
%         Vout = Vout';
%         Fout = Fout';





%
%                     t = unique(FaceIDsTemp);
%                     VisualizeTriangulation(F, V);
%                     hold on
%                     for jj = 1: length(t)
%                         fcur = t(jj);
%                         for jjj = 1:3
%                             vc = V(:, F(jjj, fcur));
%                             nc = N(:, F(jjj, fcur));                       
%                             quiver3(vc(1), vc(2), vc(3), nc(1), nc(2), nc(3));
%                             hold on
%                         end
%                         a = 2;
%                     end
%                     axis equal
%                     a = 2;
%                     for jj = 1: size(pointsTemp, 2)
%                         pc = pointsTemp(:, jj);
%                         nc = NormalsTemp(:, jj);
%                         quiver3(pc(1), pc(2), pc(3), nc(1), nc(2), nc(3));
%                         hold on
%                     end
%                     axis equal
%                     a = 2;
                    

    %                 %sort according to a distance
    %                 [dists, ix] = sort(dists, 'ascend');
    %                 pointsTemp  = pointsTemp(:, ix);
    %                 NormalsTemp = NormalsTemp(:, ix);

                    %% TOO SLOW

    %                 for j = 1:size(pointsTemp, 2)
    %                     err = err + abs(angleVec(Ncentre, NormalsTemp(:,j)));
    %                 end                 
    %                 err = err / length(ids);
    



% function pIdx = FindPointIndex(faceID, pointID)
%     pIdx = find(pointIndexing(:, 1) == faceID & pointIndexing(:, 2) == pointID);
% end


%           try
%             if stopInference(i)
%                 continue;
%             end

%             a = VAll{i};
%             b = NAll{i};
%             b = b*0.02;
%             figure();
%             VTemp = [V(:, F(1,i)), V(:, F(2,i)), V(:, F(3,i))];
%             
%             NTemp = [N(:, F(1,i)), N(:, F(2,i)), N(:, F(3,i))];
%             NTemp = NTemp*0.1;
%             
%             scatter3(a(:,1), a(:,2), a(:,3));
%             axis equal
%             hold on
%             scatter3(VTemp(1,:), VTemp(2, :), VTemp(3,:), 'red');
%             hold on
%             quiver3(a(:,1), a(:,2), a(:,3), b(:,1), b(:,2), b(:,3), 'color', 'blue');           
%             hold on
%             quiver3(VTemp(1,:), VTemp(2, :), VTemp(3,:), NTemp(1,:), NTemp(2, :), NTemp(3,:), 'color', 'blue');


%             perform_point_picking(V, F);
%         
%             trisurf(F(:, iidd)',V(1,:),V(2,:),V(3,:), 'FaceColor', [0.1, 0.1, 0.1], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 1.0);
%             light('Position',[-1.0,-1.0,100.0],'Style','infinite');
%             axis equal;
%             lighting phong;
        %     hold on
        %     plotFrame(W, vecLen, vectColors, origin); 
        
  
        
%                 dists2 = sqrt(sum((bsxfun(@minus, points, Vcentre)).^2, 1))';
%                 points = points';
%                 dists = sqrt(sum(([points(:,1)-Vcentre(1),points(:,2)-Vcentre(2),points(:,3)-Vcentre(3)]).^2, 2));
%                 dists4 = sqrt(sum((bsxfun(@minus, points, Vcentre')).^2, 2));
%                 points = points';



% %         rr = load('Temp/faceIDsTrial.mat');
% %         facesId = rr.facesId;
% %         
% %         aa = load('Temp/structCS.mat');
% %         structCS = aa.structCS;
% %         
% %         color = structCS(:, facesId);
% %         colorsLen = sqrt(color(1,:).^2 + color(2,:).^2 + color(3,:).^2);
% %         color(1,:) = color(1,:)./colorsLen;
% %         color(2,:) = color(2,:)./colorsLen;
% %         color(3,:) = color(3,:)./colorsLen;
% %         
% %         trisurf(F(:, facesId)',V(1, :),V(2, :),V(3, :), 'FaceColor', 'flat', 'FaceVertexCData', color', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
% %         hold on
% %         
% %         V_all = [];
% % 
% %         for kk = 1:length(facesId)
% % %             ids = partIDs(:, 2) == facesId(kk) & partIDs(:, 1) == 1 & partIDs(:, 3) >= 2;
% % %             V_all = [V_all; Vout(ids, :)];
% %             ids = areCentral{facesId(kk)} >= 3;
% %             V_all = [V_all; VAll{facesId(kk)}(ids, :)];
% %         end
% %         
% %         scatter3(V_all(:, 1), V_all(:, 2), V_all(:,3));




