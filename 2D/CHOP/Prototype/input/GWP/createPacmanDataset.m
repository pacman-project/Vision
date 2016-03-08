% Creates a subset of the pacman shape dataset for desired number of images
% and classes.
function [ output_args ] = createPacmanDataset( classCount, imageCount, trainingModelCount, classOrder)
     naturalPose = [3 3];
     numberOfModels = 20;
     distance = '600';
     classList = { 'Bottle', 'Bowl', 'Box', 'Can', 'Cup', 'FryingPan', 'Mug', 'TeaCup', 'Plate', 'Jug', 'TeaPot', 'ChoppingBoard', 'Shaker', 'Tray', 'Fork', 'Knife', 'Scissors', 'Spatula', 'Spoon', 'Vase'};
     modelsToAvoidTraining = [203, 213, 215, 217, 218, 219, 220];
     load([pwd '/400Pacmans/categoryNames.mat']);
     
     % Shuffle class list if needed.
     if strcmp(classOrder, 'shuffle')
          classList = classList(randperm(numel(classList)));
     end
     
     % Generate sample views, and perfect view.
     poses= allcomb(1:32, 1:8);
     numberOfPoses = size(poses,1);
     naturalPoseIdx = find(ismember(poses, naturalPose, 'rows'));
     
     % Obtain the poses we will work with.
     if imageCount == numberOfPoses
          poseIdx = (1:numberOfPoses)';
     elseif imageCount == 1
          poseIdx = naturalPoseIdx;
     else
          poseIdx = randperm(numberOfPoses, imageCount);
     end
     numberOfPoses = numel(poseIdx);
     
     % Go through every category, and create dataset.
     for classItr = 1:classCount
          trueModelIdx = find(cellfun(@(x) ~isempty(x), strfind(categoryNames, classList{classItr})));
          
          shuffledModels = randperm(numberOfModels);
          shuffledTrueModelIdx = trueModelIdx(shuffledModels);
          validTrainingModels = ~ismember(shuffledTrueModelIdx, modelsToAvoidTraining);
          if nnz(validTrainingModels) < trainingModelCount
               curTrainingModelCount = nnz(validTrainingModels);
          else
               curTrainingModelCount = trainingModelCount;
          end
          validShuffledModels = shuffledModels(validTrainingModels);
          trainingModels = validShuffledModels(1:curTrainingModelCount);
          
          % Go through the models.
          for modelItr = 1:numberOfModels
               if ismember(modelItr, trainingModels)
                    splitName = 'vocab';
                    modelId = trueModelIdx(modelItr);
               else
                    splitName = 'test';
                    modelId = trueModelIdx(modelItr);
               end
               
               inputFolder = [pwd '/400Pacmans/images' distance '/D_' num2str(modelId) '/'];
               modelFolder = [pwd '/' splitName '/' classList{classItr} '/D_' num2str(modelId) '/'];
               mkdir(modelFolder);
               
               for poseItr = 1:numberOfPoses
                    fileName = sprintf('%d_P%d_IP%d.png', modelId, poses(poseIdx(poseItr),1), poses(poseIdx(poseItr),2));
                    img = imread([inputFolder fileName]);
                    imwrite(img, [modelFolder fileName]);
               end
          end          
     end
end

