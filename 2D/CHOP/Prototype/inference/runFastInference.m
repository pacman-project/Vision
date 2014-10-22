%> Name: singleTestImages
%>
%> Description: Process each image, run discovery with vocabulary at each
%> level, and find the class each image belongs to. 
%>
%> @param testFileImages The test image names to work on.
%> @param vocabulary
%> @param redundantVocabulary
%> @param modes
%> @param options Program options.
%>
%> @retval totalInferenceTime The amount of time spent in inference.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.12.2013
function [totalInferenceTime] = runFastInference(testFileName, vocabulary, redundantVocabulary, modes, options)
    totalInferenceTime = 0;
    inferenceTimeLimit = options.inferenceTimeLimit;
    %% Get the first level nodes.
    % First, downsample the image if it is too big.
    img = imread(testFileName);
    [~, fileName, ~] = fileparts(testFileName);
    if options.debug
        display(['Processing ' fileName '.']);
    end
    % resize image if necessary.
    if max(size(img)) > options.maxImageDim
       img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
    end
    imwrite(img, [options.processedFolder '/' fileName '.png']);

    %% Form the first level nodes.
    [cellNodes, ~] = getNodes(img, [], options);
%    imwrite(smoothedImg, [options.smoothedFolder '/' fileName '.png']);
    if isempty(cellNodes)
        return;
    end
    % Save smoothed image.
    % Assign nodes their image ids.
    nodes = cell2mat(cellNodes);
    labelIds = nodes(:,1);
    
    %% Starting with the leaf nodes, we make use of indexing to speed up the process.
    % Instead of exploring all possibilities in the search space, the
    % algorithm infers high level structures on the fly, and back-tracks if
    % it cannot find anything. We put a time limit on the inference, so
    % when the time is up, it returns the highest level nodes it could
    % find. 
    startTime = tic;
    currentModes = modes{1};
    currentLevel = vocabulary{1};
    nextLevel = vocabulary{2};
    neighborhood = fix(options.edgeRadius)^2;
    inferenceTimeLimit = 10000;
    while toc(startTime) < inferenceTimeLimit
        numberOfNodes = size(nodes,1);
        childrenArr = cell(numberOfNodes,1);
        for nodeItr = 1:numberOfNodes
            centerLabelId = nodes(nodeItr,1);
            parents = currentLevel(nodes(nodeItr,1)).parents;
            nodeCoords = nodes(nodeItr,2:3);
            for parentItr = parents
                parentNode = nextLevel(parentItr);
                nodeAdjInfo = parentNode.adjInfo;
                if parentNode.children(1) ~= nodes(nodeItr,1) || nnz(~ismember(parentNode.children, labelIds)) 
                    continue;
                end
                
                % TODO: This part assumes that each composition has at
                % least 2 children! Fix if needed.
                numberOfEdges = numel(parentNode.children);
                validNode = true;
                for secNodeItr = 2:numberOfEdges
                    secNodeLabelId = parentNode.children(secNodeItr);
                    if centerLabelId == secNodeLabelId
                        validNodeIdx = ones(numberOfNodes,1)>0;
                        validNodeIdx(nodeItr) = 0;
                        validNodeIdx = validNodeIdx & labelIds == secNodeLabelId;
                    else
                        validNodeIdx = labelIds == secNodeLabelId;
                    end
                    secNodeCoords = nodes(validNodeIdx, 2:3);
                    secNodeCoords = secNodeCoords - repmat(nodeCoords, size(secNodeCoords,1),1);
                    secNodeCoords = secNodeCoords(sum(secNodeCoords.^2,2) <= neighborhood,:);
                    
                    if isempty(secNodeCoords) 
                       validNode = false;
                       break;
                    end
                    % Check for 
                    modeIdx = currentModes(:,1) == nodes(nodeItr,1) & currentModes(:,2) == secNodeLabelId;
                    testModes = currentModes(modeIdx,3:4);
                    modeIdx = find(modeIdx);
                    modeStartIdx = min(modeIdx)-1;
                    modeIdx = modeIdx - modeStartIdx;
                    expectedMode = nodeAdjInfo((secNodeItr-1),3) - modeStartIdx;
                    validMode = false;
                    for secNodeInstItr = 1:size(secNodeCoords,1)
                        [~, minIdx] = min(sum((repmat(secNodeCoords(secNodeInstItr,:), numel(modeIdx),1) - testModes).^2,2));
                        if expectedMode == minIdx
                            validMode = true;
                            break;
                        end
                    end
                    if ~validMode
                       validNode = false;
                       break;
                    end
                end
               if validNode 
                   1
               end
            end
        end
    end
    
    
    %% Process mainGraph to export realizations in the desired format for inte2D/3D integration.
    exportArr = exportRealizations(mainGraph); %#ok<NASGU>
    if exist([options.testInferenceFolder '/' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr', '-append');
    else
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr');
    end
end