%> Name: prepareData
%>
%> Description: Prepares datasets so that vocabulary learning, training and
%> testing images are separated from each other into different folders.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param fileList File list including all images in the dataset.
%> @param imageExtension The extension of the files to work on. Examples
%> include '.jpg', '.png', '_crop.png'...
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 05.12.2013 Various parameter additions, 'mode' changes
function [ ] = prepareData()
    currentPath = pwd;
    
    % Parameters to divide datasets into different sets.
    setArr = [1 2 3 4 5 10 15 20 30];
    numberOfInferenceImages = 10;
    setOffset = 0;
    
    %% Prepare sets for experiments on different classes.
    datasetFolder = [currentPath '/input/mpeg7/vocab'];
    
    fileList = fuf([datasetFolder '/*.gif'], 1, 'detail');
    vocabImagesFolder = [currentPath '/input/ECCV2014/class'];
    if exist(vocabImagesFolder, 'dir')
        rmdir(vocabImagesFolder, 's');
    end
    mkdir(vocabImagesFolder);
    for setItr = 1:numel(setArr)
        if exist([vocabImagesFolder '/' num2str(setArr(setItr)) '/vocab'], 'dir')
           rmdir([vocabImagesFolder '/' num2str(setArr(setItr)) '/vocab'], 's');
        end
        mkdir([vocabImagesFolder '/' num2str(setArr(setItr)) '/vocab']);
    end
    if exist([vocabImagesFolder '/test'], 'dir')
       rmdir([vocabImagesFolder '/test'], 's');
    end
    mkdir([vocabImagesFolder '/test']);
    
    % Mpeg-specific parameters.
    classArr = cell(70,1);
    imagePerCategory = 5;
    categoryDelim = '-';
    lastClass = '';
    classItr = 1;
    
    %% Pick samples from each class and form test cases.
    for fileItr = 1:numel(fileList)
       [~, fileName, ~] = fileparts(fileList{fileItr});
       imageClassIdx = strfind(fileName, categoryDelim);
       imageClass = fileName(1,1:(imageClassIdx(1)-1));

       if ~strcmp(imageClass, lastClass)
           lastClass = imageClass;
           classArr(classItr) = {imageClass};
           classItr = classItr + 1;
       end
    end
    testImageCount = numberOfInferenceImages;
    for classItr = 1:max(setArr)
       %% Assign some random images from each category as the vocabulary learning dataset.
       imageIdx = randsample(20, imagePerCategory+1);
       for fileItr = 1:imagePerCategory
           imagePath = [datasetFolder '/' classArr{classItr} '-0' num2str(imageIdx(fileItr)) '.gif'];
           if ~exist(imagePath, 'file')
                imagePath = [datasetFolder '/' classArr{classItr} '-' num2str(imageIdx(fileItr)) '.gif'];
           end
           validSets = setArr(setArr>=classItr);
           for validSetItr = 1:numel(validSets)
               newImagePath = [vocabImagesFolder '/' num2str(validSets(validSetItr + setOffset)) '/vocab/' num2str(classItr) '_' num2str(fileItr) '.png'];
               copyfile(imagePath, newImagePath);
           end
       end
       if testImageCount>0
           fileItr = imagePerCategory+1;
           imagePath = [datasetFolder '/' classArr{classItr} '-0' num2str(imageIdx(fileItr)) '.gif'];
           if ~exist(imagePath, 'file')
                imagePath = [datasetFolder '/' classArr{classItr} '-' num2str(imageIdx(fileItr)) '.gif'];
           end
           newImagePath = [vocabImagesFolder '/test/' num2str(testImageCount) '.png'];
           copyfile(imagePath, newImagePath);
           testImageCount = testImageCount-1;
       end
    end
    
    %% Move on to multi-view data generation.
%    setOffset = setOffset + numel(setArr);
    datasetFolder = [currentPath '/datasets/washington/coffee_mug/coffee_mug_1'];
    vocabImagesFolder = [currentPath '/input/ECCV2014/mview'];
    
    fileList = cell(max(setArr),1);
    for fileItr = 1:numel(fileList)
        fileList(fileItr) = {[datasetFolder '/coffee_mug_1_1_' num2str(round(6*fileItr)) '_crop.png']};
    end
    % Create training folders.
    for setItr = 1:numel(setArr)
        if exist([currentPath '/input/ECCV2014/mview/' num2str(setArr(setItr)) '/vocab'], 'dir') 
            rmdir([currentPath '/input/ECCV2014/mview/' num2str(setArr(setItr)) '/vocab'], 's');
        end
        mkdir([currentPath '/input/ECCV2014/mview/' num2str(setArr(setItr)) '/vocab']);
    end
    % Create test folder.
    if exist([currentPath '/input/ECCV2014/mview/test'], 'dir') 
        rmdir([currentPath '/input/ECCV2014/mview/test'], 's');
    end
    mkdir([currentPath '/input/ECCV2014/mview/test']);
    
    % Print training images.
    for fileItr = 1:numel(fileList)
        validSets = setArr(setArr>=fileItr);
        for validSetItr = 1:numel(validSets)
            newImagePath = [vocabImagesFolder '/' num2str(validSets(validSetItr + setOffset)) '/vocab/' num2str(fileItr) '.png'];
            copyfile(fileList{fileItr}, newImagePath);
        end
    end
    
    % Print test images.
    fileList = cell(numberOfInferenceImages,1);
    for fileItr = 1:numel(fileList)
        fileList(fileItr) = {[datasetFolder '/coffee_mug_1_1_' num2str(round(6*fileItr+3)) '_crop.png']};
    end
    for fileItr=1:numel(fileList)
        newImagePath = [vocabImagesFolder '/test/' num2str(fileItr) '.png'];
        copyfile(fileList{fileItr}, newImagePath);
    end
    
    %% Finally, single-category data generation.
    vocabImagesFolder = [currentPath '/input/ECCV2014/object'];
    datasetFolder = [currentPath '/input/apple/vocab'];
    fileList = fuf([datasetFolder '/*.mask.0.png'], 1, 'detail');
    
    
    for setItr = 1:numel(setArr)
        if exist([currentPath '/input/ECCV2014/object/' num2str(setArr(setItr)) '/vocab'], 'dir') 
            rmdir([currentPath '/input/ECCV2014/object/' num2str(setArr(setItr)) '/vocab'], 's');
        end
        mkdir([currentPath '/input/ECCV2014/object/' num2str(setArr(setItr)) '/vocab']);
    end
    
    for fileItr = 1:max(setArr)
        validSets = setArr(setArr>=fileItr);
        for validSetItr = 1:numel(validSets)
            newImagePath = [vocabImagesFolder '/' num2str(validSets(validSetItr + setOffset)) '/vocab/' num2str(fileItr) '.png'];
            copyfile(fileList{fileItr}, newImagePath);
        end
    end
    
    if exist([currentPath '/input/ECCV2014/object/test'], 'dir') 
        rmdir([currentPath '/input/ECCV2014/object/test'], 's');
    end
    mkdir([currentPath '/input/ECCV2014/object/test']);
    
    for fileItr = (max(setArr)+1):(max(setArr)+numberOfInferenceImages)
        newImagePath = [vocabImagesFolder '/test/' num2str(fileItr-max(setArr)) '.png'];
        copyfile(fileList{fileItr}, newImagePath);
    end
end