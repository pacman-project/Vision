%> Name: findEdgeDistanceMatrix
%>
%> Description: This function finds the pair-wise edge transformation
%> matrix. Based on the edge quantization parameter, the receptive field is 
%> divided into non-overlapping bins, and a pairwise transformation cost of
%> moving a node from one bin to another is calculated for every such pair. 
%> This pre-calculation speeds up the matching system.
%>
%> @param edgeQuantize The receptive field is divided into edgeQuantize x
%> edgeQuantize bins, uniformly distributed in the square-shaped field.
%> 
%> @param edgeSimilarityAllowed If this flag is 1, real transformation cost
%> is calculated. Otherwise, every move is assigned maximum cost, except
%> for keeping the node in the same bin (transformation with same labels).
%>
%> @retval edgeIdMatrix The edge ids of every bin.
%> @retval edgeDistanceMatrix The edge transformation matrix, which has the
%> size of nnz(edgeIdMatrix) x nnz(edgeIdMatrix). Pairwise distances
%> between all possible edge ids are calculated. 
%> @retval edgeCoords The relative coordinates in a (nnz(edgeIdMatrix) x 2)
%> array. The coordinates are relative to the center of the receptive
%> field.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 23.06.2015
function [ edgeIdMatrix, edgeDistanceMatrix, edgeCoords, logMin, logRange ] = findEdgeDistanceMatrix( edgeQuantize, distType, edgeSimilarityAllowed )
    sigma = 0.075;
    logMin = 0;
    logRange=  1;
    
    % Find a round area in which we will look for individual relations.
    edgeIdMatrix = zeros(edgeQuantize, 'int32');
    centerPointSingle = (edgeQuantize+1)/2; 
    yArr = repmat(1:edgeQuantize, edgeQuantize, 1);
    xArr = repmat((1:edgeQuantize)', 1, edgeQuantize);
    validPointMatrix = ones(edgeQuantize)>0;

    % Identify each point in the area as an individual relation.
    relations = find(validPointMatrix);
    [coordX, coordY] = ind2sub(size(validPointMatrix), relations);
    edgeCoords = [coordX, coordY] - floor(centerPointSingle);
    edgeIdMatrix(relations) = 1:numel(relations);
    
    points = [xArr(:), yArr(:)];
    
    if strcmp(distType, 'prob')
        points = (points-1) / edgeQuantize;
        endPoints = points + 1/edgeQuantize;
        topPoints = points;
        topPoints(:,1) = points(:,1) + 1/edgeQuantize;
        rightPoints = points;
        rightPoints(:,2) = points(:,2) + 1/edgeQuantize;
        
        % allocate space for distance matrix.
        distMatrix = zeros(size(points,1), 'single');
        
        % Calculate probabilities for the top left corner.
       probs = mvncdf(points, points(1,:), [sigma, sigma]);
       endProbs = mvncdf(endPoints, points(1,:), [sigma, sigma]);
       topProbs = mvncdf(topPoints, points(1,:), [sigma, sigma]);
       rightProbs = mvncdf(rightPoints, points(1,:), [sigma, sigma]);

       pointProbs = (endProbs - (topProbs + rightProbs)) + probs;

       % Create a 2D probability matrix out of the calculated probabilities.
       probMatrix = reshape(pointProbs, size(edgeIdMatrix));
       allProbMatrix = zeros(size(edgeIdMatrix)*2-1);
       allProbMatrix(edgeQuantize:end, edgeQuantize:end) = probMatrix;
       allProbMatrix(1:(edgeQuantize), edgeQuantize:end) = flipud(probMatrix);
       allProbMatrix(edgeQuantize:end, 1:(edgeQuantize)) = fliplr(probMatrix);
       allProbMatrix(1:(edgeQuantize), 1:(edgeQuantize)) = rot90(probMatrix,2);
       
       % Locate where to put the window.
       topLeftCorners = edgeCoords + centerPointSingle;
       topLeftCorners = edgeQuantize - topLeftCorners + 1;
       
       % Finally, for each location on the edge matrix, put the window in
       % the corresponding location in the probability matrix, and get
       % transition probabilities. They will be scaled for each row to add
       % up to 1. 
       for pointItr = 1:size(points,1)
          topLeftCorner = topLeftCorners(pointItr,:);
          relevantProbMatrix = allProbMatrix(topLeftCorner(1):(topLeftCorner(1) + (edgeQuantize-1)), ...
               topLeftCorner(2):(topLeftCorner(2) + (edgeQuantize-1)));
          pointProbs = relevantProbMatrix(:);
          
          % Normalize probabilities.
          pointProbs = pointProbs/sum(pointProbs);
          logProbs = single(log(pointProbs));
          distMatrix(pointItr,:) = logProbs;
       end
            
        edgeDistanceMatrix = -distMatrix(relations, relations);
        logMin = min(min(edgeDistanceMatrix));
        logRange = max(max(edgeDistanceMatrix)) - min(min(edgeDistanceMatrix));
        edgeDistanceMatrix = (edgeDistanceMatrix - logMin) / logRange;
    else
          % Find the distance matrix with respect to pairs of individual
          % relations.
          points = points(relations,:);
          distMatrix = single(pdist2(points, points, 'euclidean'));
          edgeDistanceMatrix = distMatrix / max(max(distMatrix));
    end
    
    % If edge similarities are allowed, we keep on calculating pairwise
    % edge transformation costs. Otherwise, we set them to either
    % zero-cost transformations (identical edge labels), or max-cost
    % transformations (different edge labels).
    if ~edgeSimilarityAllowed
        oneSide = size(distMatrix,1);
        edgeDistanceMatrix = inf(oneSide, oneSide, 'single');
        edgeDistanceMatrix(1:(oneSide+1):oneSide*oneSide) = 0;
    end
end

