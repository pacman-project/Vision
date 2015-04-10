%> Name: createSplits
%>
%> Description:This function creates random training splits given Caltech
%> dataset. For each split, a fixed number of test images per class is also
%> saved.
%>
%> 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
%> Ver 1.2 on 02.09.2014 Adding display commentary.
function [ ] = createSplits(datasetFolder, ext)
    numberOfRuns = 1;
%    trainCounts = [1, 3, 5, 10, 15, 20, 30];
    trainCounts = 15;
    testCount = 20; %If -1, use rest for testing.
    
    % Read names of classes.
    curFolder = pwd;
    cd(datasetFolder);
    classList = dir;
    classList = classList(3:end);
    classNames=  {classList.name};
    cd(curFolder);
    
    if ~exist([pwd '/logs'], 'dir')
        mkdir([pwd '/logs']);
    end
    
    % Read training data.
    fileNames = fuf([datasetFolder '/*' ext], 1, 'detail');
    % Assign image names into distinct classes.
    classImages = cell(numel(classNames),1);
    classImageIdx = cell(numel(classNames),1);
    for classItr = 1:numel(classImages)
        imageIdx = cellfun(@(x) ~isempty(x), strfind(fileNames, classNames{classItr}));
        classImages{classItr} = fileNames(imageIdx);
        classImageIdx{classItr} = find(imageIdx);
    end
    classImageCounts = cellfun(@(x) numel(x), classImages);
    
    % Read images into memory (so we won't have to deal with reading
    % later).
    imageArr = cell(numel(fileNames),1);
    for fileItr = 1:numel(fileNames)
        img = imread(fileNames{fileItr});
        imageArr(fileItr) = {img};
    end
    
    for runItr = 1:numberOfRuns
        for trainCount = trainCounts
            trainNames = cell(numel(classNames),1);
            testNames = cell(numel(classNames), 1);
            for classItr = 1:numel(classNames)
                % Create training and testing folders.
                classTrainFolder = [pwd '/data/split' num2str(trainCount) '/run' num2str(runItr) '/vocab/' classNames{classItr}];
                if ~exist(classTrainFolder, 'dir')
                    mkdir(classTrainFolder);
                end
                classTestFolder = [pwd '/data/split' num2str(trainCount) '/run' num2str(runItr) '/test/' classNames{classItr}];
                if ~exist(classTestFolder, 'dir')
                    mkdir(classTestFolder);
                end
                
                % Get random images for training and testing.
                numberOfImages = classImageCounts(classItr);
                randOrder = randperm(numberOfImages, numberOfImages);
                trainClassSet = randOrder(1:trainCount);
                if testCount == -1
                    classTestCount = numberOfImages-trainCount;
                else
                    classTestCount = testCount;
                end
                testClassSet = randOrder((trainCount+1):(min(numel(randOrder), classTestCount + trainCount)));
                trainNames(classItr) = {trainClassSet};
                testNames(classItr) = {testClassSet};
                
                % Write the images to the folder.
                for trainImgItr = 1:numel(trainClassSet)
                    imgOrgIdx = classImageIdx{classItr}(trainClassSet(trainImgItr));
                    [~, fileName, ~] = fileparts(fileNames{imgOrgIdx});
                    img = imageArr{imgOrgIdx};
                    imwrite(img, [classTrainFolder '/' fileName '.png']);
                end
                for testImgItr = 1:numel(testClassSet)
                    imgOrgIdx = classImageIdx{classItr}(testClassSet(testImgItr));
                    [~, fileName, ~] = fileparts(fileNames{imgOrgIdx});
                    img = imageArr{imgOrgIdx};
                    imwrite(img, [classTestFolder '/' fileName '.png']);
                end
            end
            
            % Save ids for future reference.
            save([pwd '/logs/split' num2str(trainCount) '_run' num2str(runItr) '.mat'], 'trainNames', 'testNames');
        end
    end
end

