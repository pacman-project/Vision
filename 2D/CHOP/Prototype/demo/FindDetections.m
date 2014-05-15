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
function [ finalPredictions, windowSizes ] = FindDetections( exportArr, models, feature_params, inputImgSize)
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
    windowStep = 5;
    centerDivertAllowance = 10;
    poseQuantizer = 10;
    minPredConf = 0.8; % Minimum confidence of a detection.
    detectionDivertAllowance = 20;
    
    % Update feature parameters.
    feature_params.isTesting = 1;
    
    % Find windows to process.
    allPredictions = cell(numel(windowSizes),1);
    for windowSizeItr = 1:numel(windowSizes)
        windowSize = windowSizes{windowSizeItr};
        
        dimSteps = floor((inputImgSize - windowSize)/windowStep);
        if nnz(dimSteps<0)>0 
           continue; 
        end
        
        % Find all combinations of windows.
        allCombinations = allcomb(0:dimSteps(1), 0:dimSteps(2));
        numberOfWindows = size(allCombinations,1);
        testImageSize = windowSizes{windowSizeItr};
        
        % Group realizations into separate windows
        windowExportArrs = cell(numberOfWindows,1);
        windowExportCounts = zeros(numberOfWindows,1);
        windowCenters = zeros(numberOfWindows, 2);
        for windowItr = 1:numberOfWindows
            lowXBound = allCombinations(windowItr,1)*windowStep + 1;
            lowYBound = allCombinations(windowItr,2)*windowStep + 1;
            
            windowExportArr = exportArr(exportArr(:,2) >= lowXBound & ...
                                    exportArr(:,2) < (lowXBound+windowSize(1)) & ...
                                    exportArr(:,3) >= lowYBound & ...
                                    exportArr(:,3) < (lowYBound+windowSize(2)), :);
            windowExportArr(:,2:3) = (windowExportArr(:,2:3) - ...
                    repmat([lowXBound lowYBound], size(windowExportArr,1), 1)) + 1;
                
            windowCenters(windowItr, :) = round([lowXBound lowYBound] + windowSize/2);
            windowExportArrs{windowItr} =  windowExportArr;
            windowExportCounts(windowItr) = size(windowExportArr,1);
        end
        
        % Eliminate windows based on number of realizations inside.
        validWindows = windowExportCounts >= minCounts(windowSizeItr);
        
        % Do another check based on the mean position of
        % realizations. If they do not lie close to the center, than those
        % windows can be eliminated too.
        meanPositions = cellfun(@(x) mean(x(:,2:3),1), windowExportArrs, 'UniformOutput', false);
        meanPositions = cat(1, meanPositions{:});
        validWindows = validWindows & sqrt(sum((meanPositions - repmat(windowSize/2, size(meanPositions,1), 1)).^2,2)) <= centerDivertAllowance;
        
        numberOfWindows = nnz(validWindows);
        windowExportArrs = windowExportArrs(validWindows);
        
        % Allocate space for predictions array.
        % Fill in coordinates for windows.
        % Format: [x, y, windowSizeItr, category, pose; ...]
        predictions = zeros(numberOfWindows, 6);
        predictions(:,1:2) = windowCenters(validWindows, :);
        predictions(:,3) = windowSizeItr;
        
        if numberOfWindows>0
            predictionLabels = cell(numberOfWindows,1);
            % Process each remaining window and look for detections.
            for windowItr = 1:numberOfWindows                
                % Extract features of the input image.
                [features]=featureExtractionDemo(windowExportArrs{windowItr}, 1, 1, feature_params, testImageSize);

                % Predict labels.
                [prediction] = CategoryPosePredictionDemo( features, models, 1, feature_params.integration_levels);
                category = prediction.category_prediction.labels;
                categoryConf = max(prediction.category_prediction.prob_estimates);
                pose = prediction.pose_prediction.labels;

                % Write the output.
                predictionLabels(windowItr) = {[category, mod(pose, 360), categoryConf]};
            end
            predictions(:,4:6) = cat(1, predictionLabels{:});
            allPredictions(windowSizeItr) = {predictions};
        end
    end
    allPredictions = cat(1, allPredictions{:});
    
    %% Combine overlapping predictions.
    % We combine predictions of the same type.
    detectedCategories = unique(allPredictions(:,4))';
    finalPredictions = cell(numel(detectedCategories),1);
    for categoryItr = unique(allPredictions(:,4))'
        curPredictions = allPredictions(allPredictions(:,4) == categoryItr, :);
        
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
        detImg = imdilate(detImg, strel('disk', centerDivertAllowance, 8));
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
            poses = round(poses/poseQuantizer);
            
            % Take median value of poses.
            finalPoses(detItr) = median(poses);
            
            % Learn max confidence of overlapping detections.
            confidences = confImg(detImg == detItr);
            confidences = confidences(confidences>0);
            finalConf(detItr) = mean(confidences);
        end
        
        finalPredictions(categoryItr) = {[centers(:,2), centers(:,1), ...
            repmat(curPredictions(1,3:4), size(centers,1), 1), finalPoses*poseQuantizer, finalConf]};
    end
    finalPredictions = cat(1, finalPredictions{:});
    
    %% Eliminate predictions with low probability.
    finalPredictions = finalPredictions(finalPredictions(:,6) >= minPredConf, :);
    
    %% Combine overlapping final predictions.
    if ~isempty(finalPredictions)
        [~, sortIdx] = sort(finalPredictions(:,6), 'descend');
        finalPredictions = finalPredictions(sortIdx,:);
        validDet= ones(numel(sortIdx),1)>0;
        detIdx = (1:numel(sortIdx))';
        for detItr = 1:numel(sortIdx)
            if validDet(detItr)
                overlappingDet = sqrt(sum((finalPredictions(:,1:2) - ...
                    repmat(finalPredictions(detItr,1:2), size(finalPredictions,1),1)).^2, 2)) < detectionDivertAllowance;
                validDet(overlappingDet & detIdx>detItr) = 0;
            end 
        end
        finalPredictions = finalPredictions(validDet, :);
    end
end











