% this function computes pair-wise distances between view-independant parts

function D = computePartPairWiseDistances(partsOut, layerID, pairsLeft, pairsRight, rad, isGPU)

   %  'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
   % temporary store all candidate parts as the vocabulary of the layer
   
   lim = 0.1;
   
   disp('computing pair-wise distances between parts...');
   
    lenOut = size(partsOut, 1);
    coverageOut = [];
    pairsAll = [pairsLeft; pairsRight];
    nNClusters{layerID} = lenOut;
    strOut = ['Temp/','layer', num2str(layerID), '/partSelection', num2str(layerID), '.mat'];
    save(strOut, 'partsOut', 'pairsAll', 'nNClusters', 'coverageOut');
    ReadVocabulary(layerID);
    
    Xfirst = {};
    numPoints = zeros(1, lenOut);
    vectLen = 0.01;
    computeAxis = false;
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    Norm = [0,0,1];
    Q = computeQuaternion(Norm, Norm);
    position = [0,0,0];
    
    ReadVocabulary(layerID);  % description of parts is stored to global variables

    tic
    % perform part reconsturction
    for j = 1:lenOut
        if layerID == 4 || layerID == 6
            partID = j + 1; % the first one id the planar patch
        else
            partID = j;
        end
        [partsTemp, circleRads] = VisualizePart(layerID, partID, position, Q, Norm, Xtemp, Ytemp, computeAxis, rad, vectLen, false);

%         for i = 1:size(partsTemp, 1)
%             partsTemp(i, 4:6) = partsTemp(i, 4:6)/norm(partsTemp(i, 4:6));
%         end
        
        circleIDs = zeros(size(circleRads));
        circleIDs(circleRads == 1) = 1;
        circleIDs(circleRads == 3) = 2;
        circleIDs(circleRads == 9) = 3;

        if     layerID == 3 || layerID == 4
            
            numRadialSteps = 3;
            
        elseif layerID == 5
            
            numRadialSteps = [3, 5];
            
        elseif layerID == 6
            
            numRadialSteps = [1, 3];
            
        elseif layerID == 7 || layerID == 8
            
            numRadialSteps = [0, 3, 5];
            
        elseif layerID == 9 || layerID == 10
        end
        
        
        pointsTemp = [];
        for i = 1:size(partsTemp, 1)
            NormTemp = qvrot(partsTemp(i, 4:7), [0,0,1]);
            points = circleSampling(partsTemp(i, 1:3), circleRads(i) * rad * 1.1, NormTemp, numRadialSteps(circleIDs(i)));
            pointsTemp = [pointsTemp; points];
        end
        
%         figure
%         scatter3(pointsTemp(:, 1), pointsTemp(:, 2), pointsTemp(:, 3));
%         xlim([-lim lim])
%         ylim([-lim lim])
%         zlim([-lim lim])
        
        Xfirst{j} = pointsTemp;
        numPoints(j) = size(pointsTemp, 1);
        if mod(j, 10) == 0
            disp(j);
        end
    end
    toc
    
    if isGPU
        numPointsC = cumsum(numPoints);
        Xall = vertcat(Xfirst{:});
        Xall = gpuArray(Xall);
        numPointsC = gpuArray(numPointsC);
    end
    
    % e.g  Xfirst{2} == Xall(numPointsC(1) + 1:numPointsC(2),:))
    a = 2;
    tic
    
    if ~isGPU
        tic
        D = {};
        for i = 1:lenOut
            Dtemp = zeros(1,lenOut);
            X1 = Xfirst{i};
            for j = i+1:lenOut
                dists = pdist2(X1, Xfirst{j});
                minD = min(dists, [], 1);
                Dtemp(j) = sum(minD)/length(minD);
                D{i} = Dtemp;
            end
            if mod(i, 2) == 0
                str = [num2str(i), ' out of ', num2str(lenOut)];
                disp(str);
            end
        end
        D{lenOut} = zeros(1,lenOut);
        D = vertcat(D{:});
        D = D + D';
        toc
        
    elseif isGPU % do the same on GPU
        
        idsBegin = [0, numPointsC(1:end-1)] + 1;
        idsEnd = numPointsC;
        
        tic
        D = gpuArray.zeros(lenOut, lenOut);
        
        for i = 1:lenOut
            
            Dtemp = gpuArray.zeros(1,lenOut);
            X =  Xall(idsBegin(i):idsEnd(i), :);
            Xs = sum(X.^2, 2);

            for j = i+1:lenOut
                Y = Xall(idsBegin(j):idsEnd(j), :);
                dists = sqrt(abs(bsxfun(@plus, Xs, sum(Y.^2, 2)') - 2 * (X * Y'))); % pairwise distances between points
                
%                 dists = sqrt(abs(bsxfun(@plus, Xs, Xalls(idsBegin(j):idsEnd(j))') - 2 * (X * Xall(idsBegin(j):idsEnd(j), :)'))); % pairwise distances between points
                minD = min(dists, [], 1);
                Dtemp(j) = sum(minD)/length(minD);
            end
            
            D(i, :) = Dtemp;
            if mod(i, 2) == 0
                str = [num2str(i), ' out of ', num2str(lenOut)];
                disp(str);
            end
        end
        D = gather(D);
        D(lenOut, :) = zeros(1,lenOut);
        D = D + D';
        toc
    end
    
    save('Temp/pairwiseDist.mat', 'D');
    a = 2;
end

