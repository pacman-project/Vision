% this is to perform inference of the first layer parts with FIXED size of
% the part

function PerformInference1MeshTrialFixed(list_input, lenFiles, dataSetNumber, patchRad, is_overwrite, strE)
  
    is_visualization = false;
    
    for kk = 1:lenFiles
        fullFileName = list_input{kk};        
%         fileName = 'D:\Input Data\Meshes\Aim@Shape_Selected_2.00\D00002.obj';

        [folderName, fileName] = getFolderName(fullFileName);
        outFileM = [strE,fileName(1:end-4), '.mat'];
        outFileAP = [strE, fileName(1:end-4), 'AP.mat'];

        %     if kk < 10 
        %         fileName  =  ['D0000', num2str(kk) ,'.obj'];
        %         fileNameM =  ['D0000', num2str(kk) ,'.mat'];
        %         fileNameAP = ['D0000', num2str(kk), 'AP.mat'];  % additional points
        %     elseif kk < 100 
        %         fileName  = ['D000', num2str(kk) ,'.obj'];
        %         fileNameM = ['D000', num2str(kk) ,'.mat'];
        %     elseif kk > 100
        %         fileName  = ['D00', num2str(kk) ,'.obj'];
        %         fileNameM = ['D00', num2str(kk) ,'.mat'];
        %     end

        [V, F, ~] = meshRead(fileName);

        %% make sure that size of the object is 1 at its larger dimension

        threshPlanar = 12; % degrees
        gridStep = patchRad/5;
        vecLen = 0.05;
        vectColors = eye(3);

        % visualize the mesh and some of the normals
        [V,F] = check_face_vertex(V,F);   

%         trisurf(F',V(1,:),V(2,:),V(3,:), 'FaceColor', [0.5, 0.0, 0.0], 'EdgeColor', [0.0, 0.0, 0.5], 'FaceAlpha', 1.0);
%         light('Position',[-1.0,-1.0,100.0],'Style','infinite');
%         axis equal;
%         lighting phong;
        %     hold on
        %     plotFrame(W, vecLen, vectColors, origin); 

        %     scatter3(V(1,:),V(2,:),V(3,:), 'marker', '.');
        %     axis equal;  

        lenV = size(V, 2);
        lenF = size(F, 2);

        fring = compute_face_ring(F);
        vring = compute_vertex_ring(F);
        [N,~] = compute_normal(V,F);

        % get additional points in the middle (centre of gravity) of each triangle
        %     areas = compute_Areas(V, F);
        %     q30 = quantile(areas, 0.2);
        %     q90 = quantile(areas, 0.9);

        [VAll, NAll, pointIDx, numPoints, areCentral, pointIndexing] = compute_AdditionalPointsRegular(V, F, N, gridStep);

        %     numPointsAll = sum(numPoints);
        %     points = zeros(3, numPointsAll);
        %     Normals = zeros(3, numPointsAll);
        %     startPoint = 1;
        %     
        %     for j = 1:lenF
        %         if numPoints(j) > 0
        %             endPoint = startPoint + numPoints(j) - 1;
        %             points(:,  startPoint:endPoint) = VAll{j}';
        %             Normals(:, startPoint:endPoint) = NAll{j}';
        %             startPoint = endPoint + 1;
        %         end
        %     end

        %     points = points(:, 1:endPoint);
        %     scatter3(points(1,:),points(2,:),points(3,:), 'marker', '.');
        %     axis equal; 

        %     perform_point_picking(V, F);
        % 
        %     trisurf(F(:, iidd)',V(1,:),V(2,:),V(3,:), 'FaceColor', [0.1, 0.1, 0.1], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 1.0);
        %     light('Position',[-1.0,-1.0,100.0],'Style','infinite');
        %     axis equal;
        %     lighting phong;
        %     hold on
        %     plotFrame(W, vecLen, vectColors, origin);  

        nf = 20; % can be 1:4,10

        radii = [];
        Vout = [];
        Nout = [];
        Fout = [];   % shows to which face does this part belong
        darFrames = [];
        NeiOut = {};

        parfor i = 1:lenF

            % Create a list of faces
            facesId = ExtractAdjacentFaces(i, nf, fring);

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

            for k = 1:nnz(areCentralF)

                partIDx = idsCent(k);
                Vcentre = VAll{i}(idsCent(k),:)';
                Ncentre = NAll{i}(idsCent(k),:)';

                % COMPUTE DISTANCES FROM THE CENTRAL POINT      
                dists = ComputeEuclDists(Vcentre, points, 3);
                ids = find(dists>0 & dists <= 1.5 * patchRad);
                idsN = find(dists>0 & dists <= patchRad);
                
                dists = dists(ids);
                pointsTemp = points(:, ids);
                NormalsTemp = Normals(:, ids);
                

                %sort according to a distance
                [dists, ix] = sort(dists, 'ascend');
                pointsTemp  = pointsTemp(:, ix);
                NormalsTemp = NormalsTemp(:, ix);

                rad = 0;
                counter = 0;
                is_ok = false;

                for j = 1:size(pointsTemp, 2)
                    if abs(angleVec(Ncentre, NormalsTemp(:,j))) <= threshPlanar
                        rad = dists(j);
                    else
                        counter = counter + 1;
                    end

                    if rad >= patchRad * 0.66
                        is_ok = true;
                        break;
                    end
                    if counter == 3 && ~is_ok  % no planar patch here
                        break;
                    end

                end

                if is_ok
                    radii = [radii, rad];
                    Vout = [Vout, Vcentre];
                    Nout = [Nout, Ncentre];
                    Fout = [Fout; [i, partIDx]];

                    % estimate the darboux frame here
                    [V,D] = computeDarbouxFrame(Ncentre, pointsTemp(1,:), pointsTemp(2,:), pointsTemp(3,:));
                    darFrames = [darFrames; V(:)'];
                    
                    % get a list of points that belong to this planar patch
                    FaceIDsTemp = FaceIDs(idsN);
                    pointIDsTemp = pointIDs(idsN);
                    Pind = [];
                    for nn = 1:length(FaceIDsTemp)
                        Pind = [Pind, pointIDx{FaceIDsTemp(nn)}(pointIDsTemp(nn))];
                    end
                    NeiOut{i}{partIDx} = Pind;
                end


    %             cur = cur + 1;
    %             radii(cur) = rad;
    %             Vout(:,cur) = Vcentre;
    %             Nout(:,cur) = Ncentre;
            end

            if mod(i, 100) == 0
                i
            end 
        end
    
    %     ind = find(radii > 2*patchRad/3);
    %     radii = radii(ind);
    %     Vout = Vout(:, ind);
    %     Nout = Nout(:, ind);

        lenParts = length(radii);

        % TRY TO VISUALIZE THE RECONSTRUCTION 
        V = zeros(4 * lenParts, 3);
        F = zeros(2 * lenParts, 3);
        j = 1; 

        if is_visualization
            for i = 1:length(radii)

                if radii(i) == 0
                    continue;
                end

                position = Vout(:, i);
                Norm = Nout(:, i);
                Norm = Norm/norm(Norm);

                if dot(Norm, [1;0;0]) > 10^-4
                    xtemp = [1;0;0];
                else
                    xtemp = [0;1;0];
                end

                Ytemp = cross(Norm, xtemp);
                Xtemp = cross(Norm, Ytemp);
                Ytemp = (radii(i)/2)*(Ytemp/norm(Ytemp));
                Xtemp = (radii(i)/2)*(Xtemp/norm(Xtemp));

                % establish 4 vertices and two riangles here
                V((j-1)*4 + 1, :) = (position + Xtemp + Ytemp)';  % vertex
                V((j-1)*4 + 2, :) = (position + Xtemp - Ytemp)';  % vertex  
                V((j-1)*4 + 3, :) = (position - Xtemp - Ytemp)';  % vertex  
                V((j-1)*4 + 4, :) = (position - Xtemp + Ytemp)';  % vertex  

                F((j-1)*2 + 1, :) = [(j-1)*4 + 1, (j-1)*4 + 2, (j-1)*4 + 3]; % faces
                F((j-1)*2 + 2, :) = [(j-1)*4 + 1, (j-1)*4 + 3, (j-1)*4 + 4];

                j = j+1;

                if mod(i, 1000) == 0
                    i
                end
            end

            numNormals = 1;
            rr = randperm(lenParts, numNormals);

            trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', [0.1, 0.1, 0.1], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.4);
            light('Position',[-1.0,-1.0,100.0],'Style','infinite');
            axis equal;
            lighting phong;
            hold on

            for kkk = 1:numNormals
                VD = darFrames(rr(kkk), :);
                VD = reshape(VD, [3,3]);
                Vcentre = Vout(:, rr(kkk));
                plotFrame(VD, vecLen, vectColors, Vcentre);
            end

        end  % if is_visualization

        partIDs = ones(lenParts, 3);
        partIDs(:,2:3) = Fout;

        if ~exist(strE, 'dir')
           mkdir(strE);
        end
        
        if is_overwrite || ~exist(outFileM, 'file')
            doWrite = 1;
        end
        
        if doWrite
            save(outFileM, 'Vout', 'Nout', 'darFrames', 'partIDs', 'NeiOut', 'pointIndexing');
            save(outFileAP, 'VAll', 'NAll', 'numPoints', 'areCentral');
        end
        disp(kk);
    end
end

function ed = euclDist(D1, D2)
    ed = sqrt(sum((D1 - D2).^2));
end

function an = angleVec(V1, V2) % angle of two vectors
    an = acos(dot(V1, V2) / sqrt(dot(V1, V1)*dot(V2, V2)));
    an = an * 180 / pi;
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

% function [meanV, W, Evalues, V] = PCA_model(V)
% 
%     meanV = mean(V, 2);
%     
%     V = V - repmat(meanV, [1, size(V, 2)]);
%     [W, EvalueMatrix] = eig(cov(V'));
%     Evalues = diag(EvalueMatrix);
%     [Evalues, ids] = sort(Evalues, 'descend');
%     W = W(:, ids);
%     
%     W(:,1) = W(:,1) / norm(W(:,1));
%     W(:,2) = W(:,2) / norm(W(:,2));
%     W(:,3) = W(:,3) / norm(W(:,3));
%     
%     V = W * V;
% end


%     % remove the mean variable-wise (row-wise)
%         V=V-repmat(mean(V,2),1,size(V,2));
% 
%     % calculate eigenvectors (loadings) W, and eigenvalues of the covariance matrix
%         [W, EvalueMatrix] = eig(cov(V));
%         Evalues = diag(EvalueMatrix);
% 
%     % order by largest eigenvalue
%         Evalues = Evalues(end:-1:1);
%         W = W(:,end:-1:1); W=W';  
% 
% %     % generate PCA component space (PCA scores)
%          V = W * V;



%     for i = 1:lenF
% 
%         is_planar = true;
% 
%         % set of vertices that can be fit with a planar patch
%         Vid = zeros(1,40);
%         Vcur = zeros(3,40);
%         
%         % initially we take points from 
%         Vid(1:3) = F(:, i);  
%         Vcur(:,1) = V(:,F(1, i));
%         Vcur(:,2) = V(:,F(2, i));
%         Vcur(:,3) = V(:,F(3, i));
%         
%         curCentre = Vadd(:, i);
%         normCentre = Nf(:, i);
%         planarRad = euclDist(curCentre, Vcur(:,3));
%         
%         candidates = []; % candidates to a planar patch
%         lv = 3;
% 
%         while is_planar
%             
%             if isempty(candidates)
%                 for j = 1:lv
%                     candidates = union(candidates, vring{Vid(j)}); 
%                 end
%             end
%             % circumcenter
%             
%             % sort all candidates by distances
%             dists = (V(:, candidates) - repmat(curCentre, [1, length(candidates)])).^2;
%             dists = sqrt(sum(dists,1));
%             
%             [dists, ids] = sort(dists, 'ascend');
%             candidates = candidates(ids);
%             
%             for j = 1:length(candidates) % sort candidates by distance
%                 if ismember(candidates(j), Vid)
%                     continue;
%                 else
%                     if abs(angleVec(N(:, candidates(j)), normCentre)) < threshPlanar
%                         planarRad = dists(j);
%                         lv = lv + 1;
%                         Vid(lv) = candidates(j);
%                     else
%                         is_planar = false;
%                         break;
%                     end
%                 end
%             end
%             
%             
%         end
% 
%         a = 2;
% 
%     %     V1 = V(:, curF(1));
%     %     V2 = V(:, curF(2));
%     %     V3 = V(:, curF(3));


%         for j = 1:5:numPoints*nf
%             plotFrame(Normals(:, j), vecLen, vectColors, points(:, j));
%         end
%    


%     [VCent, NCent] = compute_gravitycentre(V, F, N);

%         while length(facesId) < nf
%             if length(facesId) == 1
%                 facesId = [facesId, fring{facesId}];
%             else
%                 facesIdNew = facesId;
%                 for j = 1:length(facesId)
%                     facesIdNew = union(facesIdNew, fring{facesId(j)});
%                 end
%                 facesId = facesIdNew;
%             end
%         end



