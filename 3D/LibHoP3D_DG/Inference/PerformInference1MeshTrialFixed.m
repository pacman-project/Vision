% this is to perform inference of the first layer parts with FIXED size of
% the part

function PerformInference1MeshTrialFixed(list_input, lenFiles, patchRad, receptiveFieldRad, is_overwrite, strE)%, crossScaleStructure, curScale)
  
    is_visualization = true;
    
    kkMax = 6;
    maxNumPointDF = 100;

    scales = [0.33, 1.0, 3.0];
    
    for kk = 1:lenFiles
        
        tic
         fullFileName = list_input{kk};        
%         fullFileName = 'D:\Input Data\Meshes\cupQuadr.obj';

        [~, fileName] = getFolderName(fullFileName);
        outFileM =  [strE, '/', fileName(1:end-4), '.mat'];
        outFilePS = [strE, '/', fileName(1:end-4), 'PS.mat'];
        outFileLH = [strE, '/', fileName(1:end-4), 'LH.mat'];
        outFileNP = [strE, '/', fileName(1:end-4), 'NP.mat'];

        [V, F, ~] = meshRead(fullFileName);
%       V = perform_normal_displacement(V,F, 0.001);

        lenV = size(V, 2);
        lenF = size(F, 2);
        fringAll = ComputeFringDeep(F, kkMax, lenF);
        fring = fringAll{kkMax};
        
        [N, ~] = compute_normal(V,F);
        
        %% make sure that size of the object is 1 at its larger dimension

        threshPlanar = 10; % degrees
        gridStep = patchRad/10;
        vecLen = patchRad * 7;
        vectColors = eye(3);
        curvThresh = 1/(20*patchRad);

%         % visualize the mesh and some of the normals
%         [V,F] = check_face_vertex(V,F); 


        F1= F(1,:); F2 = F(2,:); F3 = F(3,:);
        V1 = V(:,F1); V2 = V(:,F2); V3 = V(:,F3);
%       N1 = N(:,F1); N2 = N(:,F2); N3 = N(:,F3);


        FposAll = (V1+V2+V3)/3;  % centre of each face
%         
%         N1 = N1*0.06;
%         
%         for nn = 565:570
%             
%             F(:, nn)
%             nnn = nn;
%             trisurf(F(:, nn)',V(1,:),V(2,:),V(3,:), 'FaceColor', [0.2, 0.2, 0.2], 'EdgeColor', 'black', 'FaceAlpha', 0.3);
%             light('Position',[-1.0,-1.0,100.0],'Style','infinite');
%             axis equal;
%             lighting phong;
%             axis equal;
%             hold on
% 
% %             quiver3(FposAll(1,1:nn), FposAll(2,1:nn), FposAll(3,1:nn), Nf(1,1:nn), Nf(2,1:nn), Nf(3,1:nn), 'blue');
% 
%             quiver3(V1(1,nn:nnn), V1(2,nn:nnn), V1(3,nn:nnn), N1(1,nn:nnn), N1(2,nn:nnn), N1(3,nn:nnn), 'blue');
%             a = 2;
%          end
   
%         N = N * 1.5;       
%         nn = 2740; 
%         quiver3(V(1,1:nn), V(2,1:nn), V(3,1:nn), N(1,1:nn), N(2,1:nn), N(3,1:nn), 'blue');
%         axis equal
  
        % Make a list of faces where there is no inference
%         stopInference = zeros(1, lenF);
%         for t = 1:curScale - 1
%             stopInference = stopInference + crossScaleStructure{kk}(t,:);
%         end
%         stopInference(stopInference > 1) = 1;
%         stopInference = zeros(1, lenF);

        [VAll, NAll, pointIDx, numPoints, areCentral, pointIndexing] = compute_AdditionalPointsRegular(V, F, N, gridStep); %, stopInference);
        
%         for i = 20841:20943
%             a = VAll{i};
%             figure();
%             VTemp = [V(:, F(1,i)), V(:, F(2,i)), V(:, F(3,i))];
%             scatter3(a(:,1), a(:,2), a(:,3));
%             axis equal
%             hold on
%             scatter3(VTemp(1,:), VTemp(2, :), VTemp(3,:), 'red');
%             pause(0.5);
%         end

        Vout = [];
        Nout = [];
        Fout = [];   % shows to which face does this part belong
        darFrames = [];
        likelihoods = [];
        NeiOutTemp = {};
        
        listOfFaces = zeros(1, lenF); % scale of each face
        
        for sc = 3:1
            patchRadTemp = patchRad * scales(sc);
            receptiveFieldRadTemp = receptiveFieldRad * scales(sc);
            
            parfor i = 1:lenF % 20841  

                facesId = [i, fring{i}];
                Fpos = FposAll(:,facesId)';
                Fdist = pdist2(Fpos(1,:), Fpos);

                idsF = Fdist <= 9 * receptiveFieldRadTemp;
                facesId = facesId(idsF);

                points = zeros(3, 100);
                Normals = zeros(3, 100);
                FaceIDs = zeros(1, 100);
                pointIDs = zeros(1, 100);

                % EXTRACT POINTS AND NORMALS THAT BELONG TO THESE FACES    
                startPoint = 1;
                endPoint = 0;
                for j = 1:length(facesId)
                    if numPoints(facesId(j)) > 0
                        endPoint = startPoint + numPoints(facesId(j)) - 1;
                        points(:,  startPoint:endPoint) = VAll{facesId(j)}';
                        Normals(:, startPoint:endPoint) = NAll{facesId(j)}';
                        FaceIDs(startPoint:endPoint) = facesId(j);
                        pointIDs(startPoint:endPoint) = 1:numPoints(facesId(j));
                        startPoint = endPoint + 1;
                    end
                end
                if endPoint == 0
                    continue;
                end

                points = points(:, 1:endPoint);
                Normals = Normals(:, 1:endPoint);
                FaceIDs = FaceIDs(1:endPoint);
                pointIDs = pointIDs(1:endPoint);
                areCentralF = areCentral{i};
                idsCent = find(areCentralF == 1);
                numM = nnz(areCentralF);

                for k = 1:numM

                    pointID = idsCent(k);

    %                 if partIDx == 4
    %                     a = 2;
    %                 end
                    Vcentre = VAll{i}(idsCent(k),:)';
                    Ncentre = NAll{i}(idsCent(k),:)';

                    % COMPUTE DISTANCES FROM THE CENTRAL POINT      
    %               dists = ComputeEuclDists(Vcentre, points, 3);
                    dists = sqrt(sum((points - repmat(Vcentre, 1, size(points,2))).^2, 1))';

    %                 dists2 = sqrt(sum((bsxfun(@minus, points, Vcentre)).^2, 1))';
    %                 points = points';
    %                 dists = sqrt(sum(([points(:,1)-Vcentre(1),points(:,2)-Vcentre(2),points(:,3)-Vcentre(3)]).^2, 2));
    %                 dists4 = sqrt(sum((bsxfun(@minus, points, Vcentre')).^2, 2));
    %                 points = points';

                    ids = find(dists>0 & dists <= patchRadTemp);
                    idsDF  = find(dists>0.8*patchRadTemp & dists <= 3 * receptiveFieldRadTemp);
                    idsDF2 = find(dists>2 * patchRadTemp & dists <= 9 * receptiveFieldRadTemp);

                    if length(ids) < 4 || length(idsDF) < 6
                        continue;
                    end

                    % this is done in order to speed up
                    if length(idsDF) > maxNumPointDF
                        idsTemp = randperm(length(idsDF), maxNumPointDF);
                        idsDF = idsDF(idsTemp);
                    end
                    if length(idsDF2) > maxNumPointDF2
                        idsTemp = randperm(length(idsDF2), maxNumPointDF2);
                        idsDF2 = idsDF2(idsTemp);
                    end

    %                 if length(ids) > maxNumPointsPatch
    %                     idsTemp = randperm(length(ids), maxNumPointsPatch);
    %                     ids = ids(idsTemp);
    %                 end

    %                 dists = dists(ids);
                    pointsTemp = points(:, ids);
                    NormalsTemp = Normals(:, ids);

    %                 %sort according to a distance
    %                 [dists, ix] = sort(dists, 'ascend');
    %                 pointsTemp  = pointsTemp(:, ix);
    %                 NormalsTemp = NormalsTemp(:, ix);

                    %% TOO SLOW

    %                 for j = 1:size(pointsTemp, 2)
    %                     err = err + abs(angleVec(Ncentre, NormalsTemp(:,j)));
    %                 end                 
    %                 err = err / length(ids);

                    %% TRY TO SPEED UP

                    nc = repmat(Ncentre, [1, size(pointsTemp, 2)]);
                    an = abs( acos(sum(nc.*NormalsTemp,1))  * 180 / pi); % dot(nc, NormalsTemp, 1)
                    err = sum(an) / size(pointsTemp, 2);


                    if err < threshPlanar

                        likelihood = max(0, 1 - err/threshPlanar);
                        Vout = [Vout, Vcentre];
                        Nout = [Nout, Ncentre];
                        Fout = [Fout, [i; pointID]];
                        likelihoods = [likelihoods; likelihood];

                        pointsDFTemp = points(:, idsDF);
                        normalsDFTemp = Normals(:, idsDF);
                        pointsDFTemp = pointsDFTemp - repmat(Vcentre, [1, length(idsDF)]);
                        [VV, ~] = computeDarbouxFrame_V_N(Ncentre, pointsDFTemp, normalsDFTemp, curvThresh);


                        darFrames =  [darFrames,  VV(:)];

    %                     plotFrame(VV, vecLen, vectColors, Vcentre);
    %                     hold on

                        % get a list of points that belong to this planar patch
                        FaceIDsTemp = FaceIDs(ids);
                        pointIDsTemp = pointIDs(ids);
                        Pind = pointIDx{i}(pointID);
                        for nn = 1:length(FaceIDsTemp)
                            Pind = [Pind, pointIDx{FaceIDsTemp(nn)}(pointIDsTemp(nn))];
                        end
                        NeiOutTemp{i}{pointID} = Pind;
                        %NeiOut{pointIDxCur} = Pind;
                    end

                end

                if mod(i, 1000) == 0
                    str = [num2str(i), ' out of ', num2str(lenF)];
                    disp(str);
                end 
            end
            
        end % end scale loop
        
        NeiOut = {};
        lenNO = length(NeiOutTemp);
        for i = 1:size(pointIndexing, 1)
            temp = pointIndexing(i, :);
            if temp(1) < lenNO && temp(2) < length(NeiOutTemp{temp(1)})
                NeiOut{i} = NeiOutTemp{temp(1)}{temp(2)};
            else
                NeiOut{i} = [];
            end
                
        end
        
        darFrames = darFrames';
        darFrames2 = darFrames2';
        Fout = Fout';

        % VISUALIZE THE RECONSTRUCTION 
        lenParts = length(likelihoods);
        radii = patchRad * ones(size(likelihoods));

        if is_visualization
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
            
            numNormals = 200;
            rr = randperm(lenParts, numNormals);

            for kkkk = 1:numNormals
                VD = darFrames2(rr(kkkk), :);
                VD = reshape(VD, [3,3]);
                Vcentre = Vout(:, rr(kkkk));
                plotFrame(VD, vecLen, vectColors, Vcentre);
            end

        end  % if is_visualization

        partIDs = ones(lenParts, 3);
        
        if ~isempty(partIDs)
            partIDs(:,2:3) = Fout;
        end

        if ~exist(strE, 'dir')
           mkdir(strE);
        end
        
        if is_overwrite || ~exist(outFileM, 'file')
            doWrite = 1;
        end
        if doWrite
            save(outFileM,  'Vout', 'Nout', 'darFrames', 'darFrames2', 'partIDs', 'pointIDx');  % needed for statistics collection
            save(outFilePS, 'NeiOut');                     % needed for part selection
            save(outFileLH, 'likelihoods', 'partIDs');                      % needed for  reviseFaces
            save(outFileNP, 'numPoints');                                   % needed for  reviseFaces
        end
        disp(kk);
        toc
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

function pIdx = FindPointIndex(faceID, pointID)
    pIdx = find(pointIndexing(:, 1) == faceID & pointIndexing(:, 2) == pointID);
end


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

