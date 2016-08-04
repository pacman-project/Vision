classdef PartVisualizer < handle
    %PARTVISUALIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Properties relevant to current selection.
        layerID@uint8 %layer id of the part.
        partID@int32 % part id.
        partDescription@uint8; % The modal part description image.
        minActivation@single; % Min instance activation (log probability).
        maxActivation@single; % Max instance activation (log probability).
        activationThreshold@single; % Current activation threshold.
        shareability@int32; % Part shareability.
        MDL@single; % MDL value of the part.
        numberOfInstances@int32; % Number of instances.
        keyInfo@char; % Key part info string.
        
        % Instance related info.
        instanceOffset@int32; % The current instance offset to be shown on screen.
        instanceImgs@cell; % Instance images in cell.
        instances@int32; % The instance list.
        instanceActivations@single; % Instance activations.
        instanceClasses@int32; % Instance classes.
        instanceOrder@int32; % Order of showing on screen.
        instanceLeafNodes@cell; % Instance leaf nodes.
        classMembershipArr@single; % The array for class membership contributions.
        
        % Vocabulary related stuff.
        classNames@cell; % The category name array.
        vocabulary@cell; % Vocabulary of parts.
        vocabularyDistributions@cell; % Spatial distributions of parts (Future).
        exportArr@int32; % Exported instances.
        activationArr@single; % Activations of exported instances.
        datasetName@char; % Dataset name.
        classIds@int32; % Class ids of images.
        numberOfClasses@double; % Number of classes.
        numberOfLayers@uint8; % Number of layers.
        leafNodeArr@cell; % Array containing leaf node associations.
        leafNodes@int32; % Leaf node array.
        options; % Program options.
        visFilters; % Filters for visualization.
        trainingFileNames@cell; % Name of the original files.
    end
    
    methods
        % Class constructor.
        function obj = PartVisualizer(datasetName)
            obj.layerID = uint8(1);
            obj.partID = int32(1);
            obj.instanceOffset = int32(1);
            obj.minActivation = single(-10);
            obj.maxActivation = single(0);
            obj.activationThreshold = single(0);
            
            % Read all relevant vocabulary info.
            obj = obj.readData(datasetName);
            
            % Fill with a dummy part.
            obj = obj.FillPartInfo(obj.layerID, obj.partID);
        end
        
        % Read relevant data structures.
        function obj = readData(obj, datasetName)
            obj.datasetName = datasetName;
            load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames', 'options'); 
            obj.vocabulary = vocabulary;
            obj.classNames = categoryNames;
            obj.options = options;
            obj.numberOfClasses = numel(categoryNames);
            load([pwd '/output/' datasetName '/distributions.mat'], 'vocabularyDistributions');
            obj.vocabularyDistributions = vocabularyDistributions;
            load([pwd '/output/' datasetName '/export.mat'], 'exportArr', 'activationArr', 'categoryArrIdx', 'trainingFileNames');
            obj.exportArr = exportArr;
            obj.activationArr = activationArr;
            obj.classIds = int32(categoryArrIdx);
            obj.trainingFileNames = trainingFileNames;
            
            % Read instance leaf node info.
            obj.numberOfLayers = uint8(numel(obj.vocabulary));
            allLeafNodes = cell(obj.numberOfLayers,1);
            load([pwd '/workspaces/' datasetName '/level1/leafNodes.mat']);
            load([pwd '/workspaces/' datasetName '/level1/firstLevelPrecisePositions.mat']);
            leafNodes(:,2:3) = int32(firstLevelPrecisePositions);
            obj.leafNodes = leafNodes;
            if ~exist([pwd '/workspaces/' datasetName '/instanceLeafNodes.mat'], 'file')
                for layerItr = 1:obj.numberOfLayers - 1
                    load([pwd '/workspaces/' datasetName '/level' num2str(layerItr) '/graphLevel.mat']);
                    leafNodes = {graphLevel.leafNodes};
                    allLeafNodes{layerItr} = leafNodes;
                end
                leafNodeArr = cat(2, allLeafNodes{:});
                save([pwd '/workspaces/' datasetName '/instanceLeafNodes.mat'], 'leafNodeArr');
            else
                load([pwd '/workspaces/' datasetName '/instanceLeafNodes.mat'], 'leafNodeArr');
            end
            obj.leafNodeArr = leafNodeArr;
            
            %% Create filters for visualization.
            visFilters = options.filters;
            visFilters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), visFilters, 'UniformOutput', false);
            for filterItr = 1:numel(visFilters)
                 visFilters{filterItr} = uint8(round(visFilters{filterItr} * 255));
            end
            visFilters = cat(3, visFilters{:});
            visFilters(visFilters<1) = 1;
            obj.visFilters = double(visFilters);
        end
        
        % Fill info for the relevant layer, part id. 
        function obj = FillPartInfo(obj, layerID, partID)
            obj.layerID = uint8(layerID);
            obj.partID = int32(partID);
            
            vocabNode = obj.vocabulary{layerID}(partID);
            if layerID == 1
                obj.partDescription = imread([pwd '/debug/' obj.datasetName '/level' num2str(layerID) '/reconstruction/' num2str(partID) '.png']);
                obj.MDL = single(0);
            else
                obj.partDescription = imread([pwd '/debug/' obj.datasetName '/level' num2str(layerID) '/modalProjection/' num2str(partID) '.png']); 
                obj.MDL = single(vocabNode.mdlScore);
            end
            obj.activationThreshold = single(0);
            
            % Fill in general data.
            instanceIdx = obj.exportArr(:,4) == layerID & obj.exportArr(:,1) == partID; %#ok<*PROPLC>
            obj.instanceActivations = obj.activationArr(instanceIdx);
            obj.minActivation = min(obj.instanceActivations);
            obj.maxActivation = max(obj.instanceActivations);
            obj.instances = obj.exportArr(instanceIdx,:);
            obj.numberOfInstances = int32(nnz(instanceIdx));
            obj.instanceClasses = obj.classIds(obj.instances(:,5));
            obj.instanceLeafNodes = obj.leafNodeArr(instanceIdx);
            
            % Calculate shareability.
            obj.shareability = int32(numel(fastsortedunique(sort(obj.instanceClasses))));
            
            % Class membership array info.
            classCounts = hist(single(obj.instanceClasses), single(1:obj.numberOfClasses));
            obj.classMembershipArr = classCounts ./ single(sum(classCounts));
            
            % Fill in for instances.
            obj.keyInfo = sprintf('MDL: %d\nShr:%d of %d\nNo. Inst.:%d', round(obj.MDL), obj.shareability, obj.numberOfClasses, obj.numberOfInstances);
            [~, instanceOrder] = sort(obj.instanceActivations, 'descend');
            obj.instanceOrder = int32(instanceOrder);
            obj.instanceOffset = int32(1);
        end
        
        function handles = ChangeThreshold(obj, newThreshold, handles)
            curThreshold = log(exp(obj.minActivation) + newThreshold * (exp(obj.maxActivation) - exp(obj.minActivation)));
            obj.activationThreshold = single(newThreshold);
            
            % Fill in general data.
            filterIdx = obj.instances(:,4) == obj.layerID & obj.instances(:,1) == obj.partID & obj.instanceActivations >= curThreshold;
            instances = obj.instances(filterIdx,:);
            obj.numberOfInstances = int32(size(instances,1));
            instanceClasses = obj.instanceClasses(filterIdx);
            if ~isempty(instances)
                obj.shareability = int32(numel(fastsortedunique(sort(instanceClasses))));
            else
                obj.shareability = int32(0);
            end
            
            % Class membership array info.
            classCounts = hist(single(instanceClasses), single(1:obj.numberOfClasses));
            obj.classMembershipArr = classCounts ./ single(sum(classCounts));
            
            % Fill in for instances.
            obj.keyInfo = sprintf(' MDL: %d\nShr:%d of %d\nNo. Inst.:%d', round(obj.MDL), obj.shareability, obj.numberOfClasses, obj.numberOfInstances);
            
            % Update GUI data.
            handles.visualizerData = obj;
            
            % Update GUI itself.
            handles.keyInfo.String = obj.keyInfo;
            handles.actSlider.Value = obj.activationThreshold;
            
            % Create bar graph for class contributions. 
            curClassNames = obj.classNames;
            curClassNames = cellfun(@(x) x(1:min(numel(x), 9)), curClassNames, 'UniformOutput', false);
            
            % Update class contribution graphs.
            axes(handles.classMembership);
            bar(1:obj.numberOfClasses, obj.classMembershipArr);
            set(gca,'xticklabel',curClassNames);
            set(gca,'xtick', 1:obj.numberOfClasses);
            set(gca,'xticklabelrotation', 90);
            xlim([0.4, obj.numberOfClasses+0.6]);
            
            %% Visualize instances.
            obj.instanceOffset = int32(1);
            handles = obj.VisualizeInstances(handles);
        end
        
        function handles = UpdateGUI(obj, handles)
            %% Set layer/part menus. 
            handles.levelMenu.String = cellfun(@(x) num2str(x), num2cell(1:numel(obj.vocabulary)), 'UniformOutput', false);
            handles.levelMenu.Value = obj.layerID;
            handles.partMenu.String = cellfun(@(x) num2str(x), num2cell(1:numel(obj.vocabulary{obj.layerID})), 'UniformOutput', false);
            handles.partMenu.Value = obj.partID;
            
            %% Show part description.
            axes(handles.description);
            imshow(obj.partDescription);
            
            %% Create bar graph for class contributions. 
            curClassNames = obj.classNames;
            curClassNames = cellfun(@(x) x(1:min(numel(x), 9)), curClassNames, 'UniformOutput', false);
            
            % Update class contribution graphs.
            axes(handles.classMembership);
            bar(1:obj.numberOfClasses, obj.classMembershipArr);
            set(gca,'xticklabel',curClassNames);
            set(gca,'xtick', 1:obj.numberOfClasses);
            set(gca,'xticklabelrotation', 90);
            xlim([0.4, obj.numberOfClasses+0.6]);
            
            %% Update key info.
            handles.keyInfo.String = obj.keyInfo;
            
            %% Set activation threhsolds.
            handles.lowActivation.String = num2str(obj.minActivation);
            handles.highActivation.String = num2str(obj.maxActivation);
            handles.actSlider.Value = obj.activationThreshold;
            
            %% Visualize instances.
            handles = obj.VisualizeInstances(handles);
        end
        
        function handles = VisualizeInstances(obj, handles)
           visualizedInstances = obj.instanceOrder(obj.instanceOffset:min(obj.numberOfInstances, (obj.instanceOffset+5)));
           curThreshold = log(exp(obj.minActivation) + obj.activationThreshold * (exp(obj.maxActivation) - exp(obj.minActivation)));
           imageSize = size(obj.partDescription,1);
           imageSize = [imageSize, imageSize];
           
           % First, fill spaces with dummy info.
           dummyImg = zeros(10,10,'uint8');
           for instanceItr = 1:6
                switch instanceItr
                    case 1
                        axes(handles.instance1);
                        imshow(dummyImg);
                        handles.Instance1Text.String = [];
                    case 2
                        axes(handles.instance2);
                        imshow(dummyImg);
                        handles.Instance2Text.String = [];
                    case 3
                        axes(handles.instance3);
                        imshow(dummyImg);                        
                        handles.Instance3Text.String = [];
                    case 4
                        axes(handles.instance4);
                        imshow(dummyImg);                        
                        handles.Instance4Text.String = [];
                    case 5
                        axes(handles.instance5);
                        imshow(dummyImg);
                        handles.Instance5Text.String = [];
                    case 6
                        axes(handles.instance6);
                        imshow(dummyImg);
                        handles.Instance6Text.String = [];
                end
           end

            %% Create bar graph for class contributions. 
            curClassNames = obj.classNames;
            curClassNames = cellfun(@(x) x(1:min(numel(x), 7)), curClassNames, 'UniformOutput', false);
           
            % Visualize instances one by one.
            for instanceItr = 1:numel(visualizedInstances)
                % Activation-based rendering.
                instanceActivation = obj.instanceActivations(visualizedInstances(instanceItr));
                if instanceActivation < curThreshold
                    continue;
                end
                imgID = obj.instances(visualizedInstances(instanceItr),5);
                nodes = obj.instanceLeafNodes{visualizedInstances(instanceItr)};
                experts = double(obj.leafNodes(nodes',1:3));

                % Center the nodes.
                midPoint = round(mean(experts(:,2:3),1));
                experts(:,2:3) = experts(:,2:3) - repmat(midPoint, size(experts,1), 1);
                experts(:,2:3) = round(experts(:,2:3) + repmat(imageSize/2, size(experts,1), 1));
                
                % Render part.
                [muImg, ~, ~] = obtainPoE(experts, [], [], imageSize, obj.visFilters, []);
                muImg = uint8(round(255*(double(muImg) / double(max(max(muImg))))));
                
                % Obtain the original image, crop it and visualize.
                orgImg = imread(obj.trainingFileNames{imgID});
                lowHalf = round(imageSize(1)/2) - 1;
                highHalf = imageSize(1) - lowHalf - 1;
                buffer1 = max(0, 1 - (midPoint-lowHalf));
                buffer2 = max(0, (midPoint + highHalf) - [size(orgImg,1), size(orgImg,2)]);
                croppedImg = zeros(imageSize(:,1), imageSize(:,2), size(orgImg,3), 'uint8');
                croppedImg((1+buffer1(1)):(end-buffer2(1)), (1+buffer1(2)):(end-buffer2(2)),:) = ...
                    orgImg((midPoint(1)+buffer1(1)-lowHalf):(midPoint(1)+highHalf-buffer2(1)), ...
                    (midPoint(2)+buffer1(2)-lowHalf):(midPoint(2)+highHalf-buffer2(2)), :);
                
                %% Combine both images.
                muImg(muImg == min(min(muImg))) = 0;
                muImg = label2rgb(muImg, 'jet', 'k');
                
                if size(croppedImg,3) == 1
                    croppedImg = cat(3, croppedImg, croppedImg, croppedImg);
                end
                muImg = uint8(round((single(muImg) + single(croppedImg)) / 2));
                
                % Create string and visualize samples.
                instanceString = ['Act:' num2str(instanceActivation, 3) ...
                            ', Class:' curClassNames{obj.instanceClasses(visualizedInstances(instanceItr))}];
                
                % Show the image in the screen.
                switch instanceItr
                    case 1
                        axes(handles.instance1);
                        imshow(muImg);
                        handles.Instance1Text.String = instanceString;
                    case 2
                        axes(handles.instance2);
                        imshow(muImg);
                        handles.Instance2Text.String = instanceString;
                    case 3
                        axes(handles.instance3);
                        imshow(muImg);                        
                        handles.Instance3Text.String = instanceString;
                    case 4
                        axes(handles.instance4);
                        imshow(muImg);                        
                        handles.Instance4Text.String = instanceString;
                    case 5
                        axes(handles.instance5);
                        imshow(muImg);
                        handles.Instance5Text.String = instanceString;
                    case 6
                        axes(handles.instance6);
                        imshow(muImg);
                        handles.Instance6Text.String = instanceString;
                end
            end
        end
    end
end