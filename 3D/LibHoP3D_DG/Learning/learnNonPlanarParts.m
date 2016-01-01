% this function learns non-planar first layer parts

% methodID == 1 : clustering is done based on values of mean and Gaussian curvatures.

function learnNonPlanarParts(list_input, lenFiles, dataSetNumber, patchRad, is_overwrite, strE, methodID)

is_visualization = false;
    
    for kk = 1:lenFiles
        
        fullFileName = list_input{kk};        

        [folderName, fileName] = getFolderName(fullFileName);
        outFileM =  [strE, '/', fileName(1:end-4), '.mat'];
        outFileAP = [strE, '/', fileName(1:end-4), 'AP.mat'];

        [V, F, ~] = meshRead(fullFileName);

        %% make sure that size of the object is 1 at its larger dimension

        gridStep = patchRad/5;
        vecLen = 0.02;
        vectColors = eye(3);
        curvThresh = 1/(20*patchRad);

        [V,F] = check_face_vertex(V,F);   
        lenV = size(V, 2);
        lenF = size(F, 2);

        fring = compute_face_ring(F);
        vring = compute_vertex_ring(F);
        [N,~] = compute_normal(V,F);

        [VAll, NAll, pointIDx, numPoints, areCentral, pointIndexing] = compute_AdditionalPointsRegular(V, F, N, gridStep);
        
        nf = 20; % can be 1:4,10

        Vout = [];
        Nout = [];
        Fout = [];   % shows to which face does this part belong
        darFrames = [];
        likelihoods = [];
        NeiOut = {};

        for i = 1:lenF

            % Create a list of adjacent faces
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

                pointID = idsCent(k);
                Vcentre = VAll{i}(idsCent(k),:)';
                Ncentre = NAll{i}(idsCent(k),:)';

                % COMPUTE DISTANCES FROM THE CENTRAL POINT      
                dists = ComputeEuclDists(Vcentre, points, 3);
                ids = find(dists>0 & dists <= patchRad);
                
%               dists = dists(ids);
                pointsTemp = points(:, ids);
                NormalsTemp = Normals(:, ids);
                
                pointsDFTemp = points(:, idsDF);    
                pointsDFTemp = pointsDFTemp - repmat(Vcentre, [1, length(idsDF)]);
                
                % estimate the Darboux frame here
                
                if methodID == 1
                    [VV, values] = computeDarbouxFrame(Ncentre, pointsDFTemp(1,:), pointsDFTemp(2,:), pointsDFTemp(3,:), curvThresh);
                end
%                 [H, K, QF1, QF2, S, V, D, is_ok, zsLocal] = paraboloidFitting(pointsDFTemp(1,:), pointsDFTemp(2,:), pointsDFTemp(3,:));
                
                
            end
        end

end

