% this is t o mesure similarities between vocabulary parts and grouping
% them according this measure

function ORTable = MakeOrNode(layerID, rad)
    
    str = [];
    if layerID == 3
        str = 'Temp/Layer3/partSelection3.mat';
        lim = 0.02;
    elseif layerID == 4
        str = 'Temp/layer4/partSelection4.mat';
        lim = 0.02;
    elseif layerID == 5
        str = 'Temp/layer5/partSelection5.mat';
        lim = 0.06;
    elseif layerID == 6
        str = 'Temp/layer6/partSelection6.mat';
        lim = 0.06;
    elseif layerID == 7
        str = 'Temp/layer7/partSelection7.mat';
        lim = 0.1;
    elseif layerID == 8
        str = 'Temp/layer8/partSelection8.mat';
        lim = 0.1;
     elseif layerID == 9
        str = 'Temp/layer9/partSelection9.mat';
        lim = 0.3;
    elseif layerID == 10
        str = 'Temp/layer10/partSelection10.mat';
        lim = 0.3;
    end
    
    ncb{3} = 4;
    
    dd = load(str); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
    is_visualization = true;
    
    lenOut = dd.nNClusters{layerID};
    partsOut = dd.partsOut;
    pairsAll = dd.pairsAll;
    nNClusters = dd.nNClusters;
    coverageOut = dd.coverageOut;
    
    vectLen = 0.01;
    circleRad = 0.005;
    computeAxis = false;
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    Norm =  [0,0,1];
    position = [0,0,0];
    
    numFeatures = 6;
    moments = zeros(lenOut, numFeatures);
    Xfirst = {};
    
    if layerID == 4 || layerID == 6 % do not consider a big palanar part
        lenOut = lenOut - 1;
    end
        
    ORTable = zeros(lenOut, 1); % each row has a format from->to, ex. 9->1
    
    % perform part reconsturction
    for j = 1:lenOut
        if layerID == 4 || layerID == 6
            partID = j + 1; % the first one id the planar patch
        else
            partID = j;
        end
        [partsTemp, circleRads] = VisualizePart(layerID, partID, position, Norm, Xtemp, Ytemp, computeAxis, circleRad, vectLen, false);
        for i = 1:size(partsTemp, 1)
            partsTemp(i, 4:6) = partsTemp(i, 4:6)/norm(partsTemp(i, 4:6));
        end
        % perform sampling on each circle
        pointsTemp = [];
        for i = 1:size(partsTemp, 1)
            points = circleSampling(partsTemp(i, 1:3), circleRads(i) * rad, partsTemp(i, 4:6), circleRads(i));
            pointsTemp = [pointsTemp; points];
        end
        
        Xfirst{j} = pointsTemp;
    end
    
    D = zeros(lenOut, lenOut);
    for i = 1:lenOut
        for j = i+1:lenOut
            X1 = Xfirst{i};
            X2 = Xfirst{j};
            dists = pdist2(X1, X2);
            minD = min(dists, [], 1);
            D1 = sum(minD);
            D(i, j) = D1;
            D(j, i) = D1;
        end
    end
    
    % after that we apply agglomerative hiearchical clustering
    
    L = tril(D);
    d = L(L>0)';
    Z = linkage(d, 'weighted'); %'average');
    
%     number_of_classes_unknown = true;
    
%     if number_of_classes_unknown
%         
%         bouldinIndex = zeros(1, lenOut);
%         bouldinIndex(1) = 100;
%         
%         for nCl = 1:lenOut
%             
%             c = cluster(Z,'maxclust', nCl);
%             sigmas = zeros(1, nCl);  % average intra-cluster distances
%             centralParts = [];
% 
%             for i = 1:nCl
%                 ids = find(c == i);
%                 centralPart = ids(1); % assuming parts are already sorted according to coverage
%                 d = D(centralPart, ids);
%                 centralParts = [centralParts, centralPart];
%                 sigmas(i) = sum(d)/length(d);
%             end
%             
%             Cs = D(centralParts, centralParts);
%             sumDB = 0;
%             for i = 1:nCl
%                 ratios = zeros(1, nCl);
%                 for j = 1:nCl
%                     if i == j
%                         continue;
%                     end
%                     ratios(j) = (sigmas(i) + sigmas(j))/(Cs(i,j));
%                 end
%                 sumDB = sumDB + max(ratios);
%             end
%             sumDB = sumDB/nCl;
%             bouldinIndex(nCl) = sumDB;
%         end
%     end
%     
%     [~, nCl_best] = min(bouldinIndex);

	nCl_best = ncb{layerID};
    
%     if number_of_classes_unknown
%         nCl = 0;
%         curPenalty = 0;
%         penalty = zeros(lenOut, 1);
%         withinClusterDistsAll = zeros(1, lenOut);
%         betweenClusterDistsAll = zeros(1, lenOut);
% 
%         while (nCl < lenOut) % && curPenalty < 10^(-5)
%             nCl = nCl + 1;
%             parts = zeros(nCl, 6);
%             c = cluster(Z,'maxclust', nCl);
%             withinClusterDists = 0;
%             betweenClusterDists = 0;
%             centralParts = [];
%             for i = 1:nCl
%                 ids = find(c == i);
%                 centralPart = ids(1); % assuming parts are already sorted according to coverage
%                 centralParts = [centralParts, centralPart];
%                 d = D(centralPart, ids);
%                 withinClusterDists = withinClusterDists + sum(d);
%             end
%             withinClusterDistsAll(nCl) = withinClusterDists;
%             betweenClusterDistsAll(nCl) = sum(sum(D(centralParts, centralParts)));
%         end
%     end
%     betweenClusterDistsAll = betweenClusterDistsAll;
%     x = 1:nCl;
%     plot(x, withinClusterDistsAll);
%     hold on
%     plot(x, betweenClusterDistsAll);
   
%     nCl_best = find(withinClusterDistsAll - betweenClusterDistsAll < 0);
%     
%     if layerID < 4
%         nCl_best = nCl_best(1) + 1;
%     elseif layerID <= 6
%         nCl_best = nCl_best(1) - 1;
%     else
%         nCl_best = nCl_best(1);
%     end
    
    c = cluster(Z,'maxclust', nCl_best);

    numInCluster = zeros(1, nCl_best);
    vectLen = 0.01;
    computeAxis = false;
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    Norm =  [0,0,1];
    position = [0,0,0];
    cl_reprAll = zeros(nCl_best, 1);
    
    for i = 1:nCl_best
        ids = find(c == i);
        ids = sort(ids, 'ascend');
        cl_repr = ids(1);
        cl_reprAll(i) = cl_repr;
        ORTable(ids) = cl_repr;
    end
    
    partsOutTemp = zeros(size(partsOut));
    addBegin = nCl_best + 1;
    ORTableTemp = zeros(size(ORTable));
    for i = 1:nCl_best
        curRepr = cl_reprAll(i);
        partsOutTemp(i, :) = partsOut(curRepr, :);
        ids = find(c == i);
        numPartsInCluster = length(ids) - 1;
        if numPartsInCluster > 0
            addEnd = addBegin + numPartsInCluster - 1;
            partsOutTemp(addBegin:addEnd, :) = partsOut(ids(2:end), :);
            ORTableTemp(addBegin:addEnd) = i;
            addBegin = addEnd + 1;
        end
    end
    ORTableTemp(1:nCl_best) = 1:nCl_best;
    
    
    ORTable = ORTableTemp;
    partsOut = partsOutTemp;
    nNClusters{layerID} = nCl_best;
    
    save(str, 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters');
 
    % visualizetion in order to check the results
    if is_visualization
        for i = 1:lenOut
            if layerID == 4 || layerID == 6
                partID = i + 1;
            else
                partID = i;
            end
            curCluster = ORTable(i);
            figure('Position',[20+60*(curCluster), 600 - 70 * numInCluster(curCluster),500,300]);
            is_ok = VisualizePart(layerID, partID, position, Norm, Xtemp, Ytemp, computeAxis, rad, vectLen);
            xlim([-lim lim])
            ylim([-lim lim])
            zlim([-lim lim])
            a = 2;
            numInCluster(curCluster) = numInCluster(curCluster) + 1;
        end
    end
   a = 2;
    
end

% function moments = computeMoment(partsTemp)
% 
%     % compute moments of the 2nd 3rd and the 4th order
%     lenOut = size(partsTemp, 1);
%     moments = zeros(1, 3);
%     centre = [0,0,0, 0,0,1];
%     
%     for i = 1:lenOut
%         curPart = partsTemp(i,:);
%         moments(1) = moments(1) + (centre(1) - curPart(1))^2 * (centre(2) - curPart(2))^2 * (centre(3) - curPart(3))^2; 
%         moments(2) = moments(2) + (centre(1) - curPart(1))^3 * (centre(2) - curPart(2))^3 * (centre(3) - curPart(3))^3; 
%         moments(3) = moments(3) + (centre(1) - curPart(1))^4 * (centre(2) - curPart(2))^4 * (centre(3) - curPart(3))^4; 
%     end
% end

%         moments(j, 1:3) = computeMoment(partsTemp(:, 1:3));
%         moments(j, 4:6) = computeMoment(partsTemp(:, 4:6));












