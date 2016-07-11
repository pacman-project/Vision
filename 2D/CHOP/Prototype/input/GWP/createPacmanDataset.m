% Creates a subset of the pacman shape dataset for desired number of images
% and classes.
function [ ] = createPacmanDataset( classCount, imageCount, trainingModelCount, testImageCount, testModelCount, classOrder)
     numberOfModels = 20;
     datasetStr = '400Pacmans';
     distance = '900';
     classList = { 'Bottle', 'Bowl', 'Box', 'Mug', 'Cup', 'FryingPan', 'Spatula', 'TeaPot', 'TeaCup', 'Plate', 'Jug', 'Can', 'ChoppingBoard', 'Shaker', 'Tray', 'Fork', 'Knife', 'Scissors', 'Spoon', 'Vase'};
     load([pwd '/' datasetStr '/categoryNames.mat']);
     
     % Shuffle class list if needed.
     if strcmp(classOrder, 'shuffle')
          classList = classList(randperm(numel(classList)));
     end
     
     % Generate sample views, and perfect view.
     numberOfPoses = 256;
     
     % Go through every category, and create dataset.
     for classItr = 1:classCount
         
          % Decide on models to train and test.
          shuffledModelIds = datasample(1:numberOfModels, numberOfModels, 'Replace', false);
          trainingModels = shuffledModelIds(1:trainingModelCount);
          testModels = shuffledModelIds((trainingModelCount+1):(trainingModelCount+testModelCount));
         
          % Go through the models.
          for modelItr = 1:numberOfModels
              
              % What's actual model id?
              trueModelIdx = find(cellfun(@(x) ~isempty(x), strfind(categoryNames, classList{classItr})));
              
              % Based on whether this specific model is train/test, set
              % paths/other stuff.
              if ismember(modelItr, trainingModels)
                   splitName = 'vocab';
                   curImageCount = imageCount;
                   modelId = trueModelIdx(modelItr);
              elseif ismember(modelItr, testModels)
                   splitName = 'test';
                   curImageCount = testImageCount;
                   modelId = trueModelIdx(modelItr);
              else
                  continue;
              end
               
              % Get poses.
              poseIdx = datasample(1:numberOfPoses, curImageCount, 'Replace', false);
              
              % Obtain image list and save new images.
              imgList = dir([pwd '/' datasetStr '/images' distance '/D_' num2str(trueModelIdx(modelItr))]);
              imgList = imgList(3:end);
              
              % Set folders!
              inputFolder = [pwd '/' datasetStr '/images' distance '/D_' num2str(modelId) '/'];
              modelFolder = [pwd '/' splitName '/' classList{classItr} '/D_' num2str(modelId) '/'];
              mkdir(modelFolder);
               
              % Write images in corresponding folders.
              for poseItr = 1:curImageCount
                   fileName = imgList(poseIdx(poseItr)).name;
                   img = imread([inputFolder fileName]);
                   imwrite(img, [modelFolder fileName]);
              end
          end          
     end
end

