% this is t o mesure similarities between vocabulary parts and grouping
% them according this measure

function [ output_args ] = MakeOrNode(layerID)

    if layerID == 3
        dd = load('Temp/Layer3/partSelection3.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.02;
    elseif layerID == 4
        dd = load('Temp/layer4/partSelection4.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.02;
    elseif layerID == 5
        dd = load('Temp/partSelection5.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.06;
    elseif layerID == 6
        dd = load('Temp/partSelection6.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.06;
    end
    lenOut = dd.nNClusters{layerID};
    
    vectLen = 0.01;
    circleRad = 0.005;
    computeAxis = false;
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    Norm =  [0,0,1];
    position = [0,0,0];
    
    
    numFeatures = 6;
    moments = zeros(lenOut, numFeatures);
    
    % perform part reconsturction
    for j = 1:lenOut
        [partsTemp, ~] = VisualizePart(layerID, j, position, Norm, Xtemp, Ytemp, computeAxis, circleRad, vectLen, false);
        for i = 1:size(partsTemp, 1)
            partsTemp(i, 4:6) = partsTemp(i, 4:6)/norm(partsTemp(i, 4:6));
        end
        moments(j, 1:3) = computeMoment(partsTemp(:, 1:3));
        moments(j, 4:6) = computeMoment(partsTemp(:, 4:6));
    end
    
    % normalize moments
    for i = 1:numFeatures
        sumM = sum(abs(moments(:, i)))/lenOut;
        moments(:, i) = moments(:, i)/sumM;
    end
    
    % apply kmeans clustering optimizing ratio
    % intra-clustre-distance/inter-cluster-dist
    ratio = zeros(1, lenOut);
    
    for i = 5 % number of clusters
        idx = kmeans(moments, i);
        clusterCentres = zeros(i, 3);
        %measure statistics of each cluster
        for j = 1:i
            ids = idx == j;
            points = moments(ids, :);
            clusterCentres(j, :) = sum(points)/sum(ids);
            % find a sample which is closest to the cluster center
            dist = pdist2(clusterCentres(j, :), points);
            idDM = find(dist == min(dist));
            idDM = idDM(1);
            clusterCentres(j, :) = points(idDM, :);
        end
    end
    a = 2;
end

function moments = computeMoment(partsTemp)

    % compute moments of the 2nd 3rd and the 4th order
    lenOut = size(partsTemp, 1);
    moments = zeros(1, 3);
    centre = [0,0,0, 0,0,1];
    
    for i = 1:lenOut
        curPart = partsTemp(i,:);
        moments(1) = moments(1) + (centre(1) - curPart(1))^2 * (centre(2) - curPart(2))^2 * (centre(3) - curPart(3))^2; 
        moments(2) = moments(2) + (centre(1) - curPart(1))^3 * (centre(2) - curPart(2))^3 * (centre(3) - curPart(3))^3; 
        moments(3) = moments(3) + (centre(1) - curPart(1))^4 * (centre(2) - curPart(2))^4 * (centre(3) - curPart(3))^4; 
    end
end












