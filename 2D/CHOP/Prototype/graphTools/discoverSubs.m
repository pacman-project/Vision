%> Name: discoverSubs
%>
%> Description: This function discovers parts and their realizations with
%> graphs defined by vocabLevel and graphLevel, as well as their printed
%> versions in graphFileName. 
%>
%> @param vocabLevel If preDefinedSearch is 1, the compositions in this
%> vocabulary level are searched in graphLevel. If 0, simply ignored.
%> @param graphLevel The current object graphs' level.
%> @param options Program options.
%> @param currentFolder Path to the workspace folder.
%> @param preDefinedSearch If true, supervised search is run. 
%> If empty, unsupervised SUBDUE runs over graphLevel.
%> @param levelItr current level id.
%>
%> @retval vocabLevel If preDefinedSearch is 1, [] is returned. If it is 0,
%> best compositions discovered are returned.
%> @retval graphLevel The graph consisting of part realizations discovered.
%> It's just a dummy graph including realization labels and children.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.01.2014
%> 'self' type search added on 05.02.2014
function [vocabLevel, graphLevel] = discoverSubs( vocabLevel, graphLevel, nodeDistanceMatrix, options, preDefinedSearch, threshold, levelItr)
    startTime = tic;
    edgeDistanceMatrix = options.edgeDistanceMatrix;
    if ~preDefinedSearch
        display(['.... Discovering compositions in level ' num2str(levelItr) '.']); 
    end
    load([options.currentFolder '/output/' options.datasetName '/export.mat'], 'categoryArrIdx');
    
    % Search for substructures.
    if preDefinedSearch
        % Inference on the test image with learned vocabularies.
        % This part is again related to combining parts. 
        % The status of this section is debatable. Please wait for updates.
        % It'll be unhid as soon as possible.
        graphLevel = inferSubs(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, threshold);
    else
        %% We're experimenting here. Instead of setting a human-set threshold for similarity, 
        % we try to limit the number of compositions at each layer.
        % Then, we search an optimal threshold value (Using sort-of a binary search).
        hiddenNodeCount = options.reconstruction.numberOfReconstructiveSubs;  
        minThr =  options.subdue.minThreshold; % Minimum threshold for elastic matching.
        maxThr = options.subdue.maxThreshold; % Max threshold for elastic part matching. 
        midThr = (minThr + maxThr) / 2;
        currentDepth = 1;
        maxDepth = options.subdue.thresholdSearchMaxDepth;
        
        display(['[SUBDUE] Starting threshold search with ' num2str(midThr) '. We are limited to ' num2str(hiddenNodeCount) ' compositions.']);
        while (currentDepth <= maxDepth)
            options.subdue.threshold = midThr;
%            [T, tmpVocabLevel, tmpGraphLevel, isCoverageOptimal] = evalc('runSubdue(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, options);');
            [tmpVocabLevel, tmpGraphLevel, isCoverageOptimal] = runSubdue(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, options);
            % Depending on which side we are on the ideal hidden node count, 
            % we make a binary search on the "good" threshold.
            if numel(tmpVocabLevel) == hiddenNodeCount && isCoverageOptimal
                % If we've found the perfect number of hidden nodes, exit.
                display(['[SUBDUE] Found perfect similarity threshold: ' num2str(midThr) '! Quitting..']);
                break;
            elseif numel(tmpVocabLevel) < hiddenNodeCount
                maxThr = midThr;
                midThr = (maxThr + minThr) / 2;
                display(['[SUBDUE] Too few generated compositions (' num2str(numel(tmpVocabLevel)) ').' ...
                    ' Lowering the threshold to ' num2str(midThr) ' and continuing to search..']);
            else
                % We have far too many hidden nodes. Make threshold
                % smaller.                              
                minThr = midThr;
                midThr = (maxThr + minThr) / 2;
                display(['[SUBDUE] Too many generated compositions (' num2str(numel(tmpVocabLevel)) ').' ...
                    ' Increasing the threshold to ' num2str(midThr) ' and continuing to search..']);
            end
            currentDepth = currentDepth + 1;
            % Open threads for parallel processing.
            if options.parallelProcessing
                s = matlabpool('size');
                if s>0
                   matlabpool close; 
                end
                matlabpool('open', options.numberOfThreads);
            end
        end
%        display(T);
        if currentDepth>maxDepth && numel(tmpVocabLevel) ~= hiddenNodeCount
            display(['[SUBDUE] Maximum depth reached. Best approximation for similarity threshold:' num2str(midThr)]);
        end
        display('[SUBDUE] Saving the vocabulary, instances and the threshold.');
        % Save vocabLevel and graphLevel.
        vocabLevel = tmpVocabLevel;
        graphLevel = tmpGraphLevel;
    end
    
    % Show time elapsed.
    display(['.... Time elapsed: ' num2str(toc(startTime)) ' secs.']);
    display(['.... Found ' ...
        num2str(numel(graphLevel)) ' instances of ' num2str(numel(vocabLevel)) ' compositions.']);
end
