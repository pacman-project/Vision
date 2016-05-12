% Creates a subset of the pacman shape dataset for desired number of images
% and classes.
function [ ] = createPacmanDataset( classCount, imageCount, trainingModelCount, classOrder)
     numberOfModels = 20;
     datasetStr = '400Pacmans';
     distance = '600';
     classList = { 'Bottle', 'Bowl', 'Box', 'Mug', 'Cup', 'FryingPan', 'Spatula', 'TeaPot', 'TeaCup', 'Plate', 'Jug', 'Can', 'ChoppingBoard', 'Shaker', 'Tray', 'Fork', 'Knife', 'Scissors', 'Spoon', 'Vase'};
     load([pwd '/' datasetStr '/categoryNames.mat']);
     
     % Shuffle class list if needed.
     if strcmp(classOrder, 'shuffle')
          classList = classList(randperm(numel(classList)));
     end
     
     % Generate sample views, and perfect view.
%     poses= allcomb(1:32, 1:8);
     poseIdx = datasample(1:256, imageCount, 'Replace', false);
     numberOfPoses = numel(poseIdx);
     
     % Go through every category, and create dataset.
     for classItr = 1:classCount
         
          % Go through the models.
          for modelItr = 1:numberOfModels
              trueModelIdx = find(cellfun(@(x) ~isempty(x), strfind(categoryNames, classList{classItr})));
              trainingModels = 1:trainingModelCount;
              
              imgList = dir([pwd '/' datasetStr '/images' distance '/D_' num2str(trueModelIdx(modelItr))]);
              imgList = imgList(3:end);
              
               if ismember(modelItr, trainingModels)
                    splitName = 'vocab';
                    modelId = trueModelIdx(modelItr);
               else
                    splitName = 'test';
                    modelId = trueModelIdx(modelItr);
               end
               
               inputFolder = [pwd '/' datasetStr '/images' distance '/D_' num2str(modelId) '/'];
               modelFolder = [pwd '/' splitName '/' classList{classItr} '/D_' num2str(modelId) '/'];
               mkdir(modelFolder);
               
               for poseItr = 1:numberOfPoses
                    fileName = imgList(poseIdx(poseItr)).name;
                    img = imread([inputFolder fileName]);
                    imwrite(img, [modelFolder fileName]);
               end
          end          
     end
end

