function [partsOutSpecial, nNClustersSpecial, partEntropySpecial] = partSelectionSpecial(is_Entropy, is_localization, partsLayerEnt, partsLayerLoc, nClusters, ...
                                         fileForVisualizationPrevLayer, clusterCurDepths, displ, ...
                                         layerID, cluster1Centres, depthStep, fieldSize)
   
%    load(statisticsLayerSieved);
%    load(statisticsLayerAggregated);  %  X, etc.
%    
%    lenCombs = size(X, 1);
   n2Clusters = nClusters ^2;
   
   downsamplingScheme = 3;
   
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    fieldCenter = halfFieldSize + 1;
   
   
   if is_Entropy
       load(partsLayerEnt);  % 'triplesOutEnt', 'numSelectedEnt'  'partEntropyOut'
   else 
       partEntropyOut = [];
       triplesOutEnt = [];
       numSelectedEnt = 0;
   end
   
   if is_localization
       cc = load(partsLayerLoc);  % 'triplesOutLoc', 'numSelectedLoc'  'partLocEntropyOut'
       triplesOutLoc = cc.triplesOutLoc;
       numSelectedLoc = cc.numSelectedLoc;
       partLocEntropyOut = cc.partLocEntropyOut;
   else
       triplesOutLoc = [];
       numSelectedLoc = 0;
       partLocEntropyOut = [];
   end

    
    numIter = [numSelectedEnt, numSelectedLoc];
    partsCand{1} = triplesOutEnt;
    partsCand{2} = triplesOutLoc;
    
    sortCriteria{1} = partEntropyOut;
    sortCriteria{2} = partLocEntropyOut;
    
    partsEntropy{1} = partEntropyOut;
    partsEntropy{2} = partLocEntropyOut;
    
    partsOutSpecial = [];
    partEntropySpecial = [];
    
    % do part selection iteratively
    for iii = 2:2
        
        str = ['part selection special ', num2str(iii)];
        disp(str);
          
        X = partsCand{iii};
        curEntropy = partsEntropy{iii};
        lenCombs = size(X, 1);
        is_selected = zeros(1, lenCombs);

        [~, ind] = sort(sortCriteria{iii}, 'ascend');
        X = X(ind, :);
        curEntropy = curEntropy(ind);
        
        % otherwise aggregate them according to their similarities
        [X_first, ~, emptyIndicator] = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer{layerID-1},  ...
                       displ, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths);  % this is a first layer descriptor

        X_first(emptyIndicator == 1) = fieldSize(3)*3*depthStep;
        distancesAll = Integral_distances(X_first, X_first, lenCombs, lenCombs, false, false);  % compute distances from all parts to 
        
        distAll = distancesAll(:);
        distAll = sort(distAll, 'ascend');
        distThresh = distAll(7*lenCombs);

        
        for i = 1:numIter(iii)
            
            if is_selected(i) 
                continue;
            end
            
            % take the first part
            partsOutSpecial = [partsOutSpecial; X(i, :)];
            partEntropySpecial = [partEntropySpecial; curEntropy(i)]; 
            dists = distancesAll(i,:);
            indss = dists<=distThresh;
            is_selected(indss) = 1;

        end
    end
    
    nNClustersSpecial = size(partsOutSpecial, 1);
    
end

