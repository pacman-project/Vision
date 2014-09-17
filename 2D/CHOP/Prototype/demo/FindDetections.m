%> Name: FindDetections
%>
%> Description: This function processes the output of test image and finds
%> its category/pose label.
%>
%> @param exportArr Exported realizations of the test image.
%> @param models Models to be used in prediction of category and pose.
%> @param feature_params Feature parameters for prediction features.
%> @param inputImgSize Size of the input image.
%>
%> @retval detections Detection array of the form:
%> [x, y, windowSizeItr, category, pose;
%> [x, y, windowSizeItr, category, pose;
%> ...]
%> @retval windowSizes Window sizes, indexed with windowSizeItr.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 14.05.2014
function [ finalPredictions, windowSizes, finalGroupMask] = FindDetections( exportArr, models, feature_params, inputImgSize)
    %% Here, we process the realizations in a window-based approach.
    % We process the coordinates using sliding windows to work in a
    % localized environment for detections. Some windows that do not have
    % enough realizations to work with are not processed for efficiency
    % reasons.
    
    % Set parameters
    trainWindowSizes = cat(1, feature_params.imageSizes{:});
    windowSizes = {[180, 220], [285, 360]};
    minCounts = zeros(numel(windowSizes), 1);
    for windowSizeItr = 1:numel(windowSizes)
        minCounts(windowSizeItr) = min(feature_params.minCounts(ismember(trainWindowSizes, windowSizes{windowSizeItr}, 'rows')));
    end
    realizationGroupingThr = 10;       % Proximity to group level 1 realizations.
    sameDetectionDivertAllowance = 10; % Maximum distance of detections 
                                       % of the same class.
    minMeanLevel = 4; % Minimum level of compositions used to validate 
                        % windows based on their mean positions.
    minPredConf = 0.0; % Minimum confidence of a detection.
    detectionDivertAllowance = 50;
    
    % Update feature parameters.
    feature_params.isTesting = 1;
    
    % We'll find the windows in a different manner. Group nearby level 1
    % nodes, and run connected component analysis to find clusters. Each
    % cluster is a possible window.
    maskImg=zeros(inputImgSize);
    % Get positions as liner indices, and mark them in a mask.
    pos = exportArr(exportArr(:,4)==1,2:3);
    ind = sub2ind(inputImgSize, pos(:,1), pos(:,2));
    maskImg(ind) = find(exportArr(:,4)==1);

    % Group the realizations based on proximity.
    groupedMask = bwlabel(imfill(imdilate(maskImg>0, strel('disk', realizationGroupingThr)), 'holes'));
        
    % Find windows to process.
    allPredictions = cell(numel(windowSizes),1);
    for windowSizeItr = 1:numel(windowSizes)
        windowSize = windowSizes{windowSizeItr};
        
        % If input image is too small, do nothing.
        if nnz(inputImgSize < windowSize)
           continue; 
        end
        
        % Allocate space for windows.
        halfSize = ceil(windowSize/2)-1;
        numberOfWindows = numel(setdiff(unique(groupedMask), 0));
        windowExportArrs = cell(numberOfWindows,1);
        windowExportCounts = zeros(numberOfWindows,1);
        windowCenters = zeros(numberOfWindows, 2);
        validWindowCtr = 0;
        
        %% Process each window separately. 
        for groupItr = 1:numberOfWindows
            % Get corresponding realizations.
            groupedMaskInd = find(groupedMask == groupItr);
            realizations = setdiff(unique(maskImg(groupedMaskInd)),0);
            
            % Find the center of the window. 
            % In small windows, average the level 1 realizations'
            % positions. In large ones, simply take the center of the
            % bounding box of the group's mask.
            if windowSizeItr == 1
                center = round(mean(exportArr(realizations, 2:3), 1));
            else
                [xInd, yInd] = ind2sub(inputImgSize, groupedMaskInd);
                center = round([(min(xInd) + max(xInd))/2, (min(yInd) + max(yInd))/2]);
            end
            
            % Check if center is not too close to the sides.
             if center(1) < 100 || center(2) < 100 || ...
                 center(1) >= (inputImgSize(1) - 100) || center(2) >= (inputImgSize(2) - 100)
                 continue;
             end
            
            % Get all realizations.
            exportArrInd = sub2ind(inputImgSize, exportArr(:,2), exportArr(:,3));
            realizations = ismember(exportArrInd, groupedMaskInd);
            windowExportArr = exportArr(realizations,:);
            prevSize = size(windowExportArr,1);
            
            % Check if all realizations lie within the mask.
            windowExportArr = windowExportArr(windowExportArr(:,2) >= (center(1) - halfSize(1))  & ...
                                 windowExportArr(:,2) <= (center(1) + halfSize(1)) & ...
                                 windowExportArr(:,3) >= (center(2) - halfSize(2)) & ...
                                 windowExportArr(:,3) <= (center(2) + halfSize(2)), :);
            newSize = size(windowExportArr,1);
            
            % If some realizations are eliminated, we need to exit.
            if newSize ~= prevSize
                continue;
            end
            
            % Increase counter, we've got a good window.
            validWindowCtr = validWindowCtr + 1;
            
            % Put the window to the center.
            windowCenters(validWindowCtr, :) = center;
            
            % Write output for this window.
            windowExportArr(:,2:3) = (windowExportArr(:,2:3) - repmat((center - halfSize), size(windowExportArr,1), 1)) + 1;
            windowExportArrs(validWindowCtr,:) = {windowExportArr};
            windowExportCounts(validWindowCtr) = size(windowExportArr,1);
        end
        
        if validWindowCtr == 0
            continue;
        end
        
        % Remove empty entries.
        windowCenters = windowCenters(1:validWindowCtr,:);
        windowExportArrs = windowExportArrs(1:validWindowCtr,:);
        windowExportCounts = windowExportCounts(1:validWindowCtr,:);

        % Eliminate windows based on number of realizations inside.
        validWindows = windowExportCounts >= (minCounts(windowSizeItr));
        
        % Eliminate windows that do not have high-level compositions.
        validWindows = validWindows & cellfun(@(x) nnz(x(:,4) >= minMeanLevel) > 0, windowExportArrs);
        
        % Save those windows which have passed all checks!
        numberOfWindows = nnz(validWindows);
        windowExportArrs = windowExportArrs(validWindows);
        
        %% Fill in coordinates and other info for each window.
        % Allocate space for predictions array.
        % Format: [x, y, windowSizeItr, category, pose, confidence; ...]
        predictions = zeros(numberOfWindows, 6);
        predictions(:,1:2) = windowCenters(validWindows, :);
        predictions(:,3) = windowSizeItr;
        
        if numberOfWindows>0
            predictionLabels = cell(numberOfWindows,1);
            % Process each remaining window and look for detections.
            for windowItr = 1:numberOfWindows   
                testImageSize = windowSizes{windowSizeItr}; 
                % Extract features of the input image.
%                [features]=featureExtractionDemo(windowExportArrs{windowItr}, 1, 1, feature_params, testImageSize);
                locations = windowExportArrs{windowItr}(:,2:3);
                if ismember(180, testImageSize)
                    testImageSize = [285, 360];
                    locations = locations + repmat([53, 70], size(locations,1), 1);
                end
                maskImg2=zeros(testImageSize);
                iddx = sub2ind(testImageSize, locations(:,1), locations(:,2));
                maskImg2(iddx)=1;
                features = HOG(maskImg2)';

                % Predict labels.
                [prediction] = CategoryPosePredictionDemo( features, models, 1);
                category = prediction.category_prediction.labels;
                categoryConf = max(prediction.category_prediction.prob_estimates);
                pose = prediction.pose_prediction.labels;

                % Write the output.
                predictionLabels(windowItr) = {[category, pose, categoryConf]};
            end
            predictions(:,4:6) = cat(1, predictionLabels{:});
            allPredictions(windowSizeItr) = {predictions};
        end
    end
    allPredictions = cat(1, allPredictions{:});
   
    %% If no predictions are found, return.
    if isempty(allPredictions)
        finalPredictions = [];
        windowSizes = [];
        finalGroupMask = [];
        return;
    end
    
    %% Combine overlapping predictions.
    % We combine predictions of the same type.
    detectedCategories = unique(allPredictions(:,4))';
    finalPredictions = cell(numel(detectedCategories),1);
    for categoryItr = unique(allPredictions(:,4))'
        curPredictions = allPredictions(allPredictions(:,4) == categoryItr, :);
        
        % Prepare helper matrices to estimate average pose and confidence.
        detImg = zeros(inputImgSize);
        poseImg = ones(inputImgSize);
        poseImg = poseImg * -1;
        confImg = zeros(inputImgSize);
        for preItr = 1:size(curPredictions,1)
            detImg(curPredictions(preItr,1), curPredictions(preItr,2)) = 1;
            poseImg(curPredictions(preItr,1), curPredictions(preItr,2)) = curPredictions(preItr,5);
            confImg(curPredictions(preItr,1), curPredictions(preItr,2)) = curPredictions(preItr,6);
        end
        
        % Dilate the image containing detection points.
        detImg = detImg > 0;
        detImg = imdilate(detImg, strel('disk', sameDetectionDivertAllowance, 8));
        detImg = bwlabel(detImg);
        
        % Find center points to report as detections.
        centers = regionprops(detImg, 'centroid');
        centers = round(cat(1, centers.Centroid));
        finalPoses = zeros(size(centers,1), 1);
        finalConf = zeros(size(centers,1), 1);
        
        % Find average poses for each separate detection
        for detItr = 1:size(centers,1)
            poses = poseImg(detImg == detItr);
            poses = poses(poses>-1);
            
            % Take median value of poses.
            finalPoses(detItr) = mod(round(median(poses)),12);
            
            % Learn mean confidence of overlapping detections.
            confidences = confImg(detImg == detItr);
            confidences = confidences(confidences>0);
            finalConf(detItr) = mean(confidences);
        end
        finalPredictions(categoryItr) = {[centers(:,2), centers(:,1), ...
            repmat(curPredictions(1,3:4), size(centers,1), 1), finalPoses, finalConf]};
    end
    finalPredictions = cat(1, finalPredictions{:});
    
    %% Eliminate predictions with low probability.
    finalPredictions = finalPredictions(finalPredictions(:,6) >= minPredConf, :);
    
    %% Combine overlapping predictions.
    if isempty(finalPredictions)
       display('Empty'); 
    else
        % Sort them based on their scores.
        [~, sortIdx] = sort(finalPredictions(:,6), 'descend');
        finalPredictions = finalPredictions(sortIdx,:);
        
        % Get overlapping windows, and let them vote for categories.
        validDetections = ones(size(finalPredictions,1),1)>0;
        finalPositions = finalPredictions(:,1:2);
        for finalItr = 1:size(finalPredictions,1)
            if validDetections(finalItr)
                distances = sqrt(sum((finalPositions - repmat(finalPositions(finalItr,:), size(finalPredictions,1), 1)).^2, 2));
                validDetections(distances <= detectionDivertAllowance & validDetections) = 0;
                validDetections(finalItr) = 1;
            end
        end
        finalPredictions = finalPredictions(validDetections,:);
    end
    
    % Return the groups that correspond to objects.
    finalGroupMask = zeros(inputImgSize)>0;
    centers = regionprops(groupedMask, 'Centroid');
    centers = cat(1, centers.Centroid);
    centers = [centers(:,2), centers(:,1)];
    for finalItr = 1:size(finalPredictions,1)
        distances = sqrt(sum((centers - repmat(finalPredictions(finalItr,1:2), size(centers,1),1)).^2,2));
        [~, groupToAdd] = min(distances);
        if groupToAdd > 0
            finalGroupMask(groupedMask == groupToAdd) = 1;
        end
    end
end