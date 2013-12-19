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
function [ ] = prepareData( datasetName, numberOfVocabImagesPerCategory)
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    % Set folder variables.
    datasetFolder = [currentPath '/datasets/' datasetName];
    
    fileList = fuf([datasetFolder '/*', '.gif'], 1, 'detail');
    vocabImageFolder = [currentPath '/input/' datasetName '/vocab'];
    if exist(vocabImageFolder, 'file')
        rmdir(vocabImageFolder, 's');
    end
    mkdir(vocabImageFolder);
    if strcmp(datasetName, 'mpeg7')
       classArr = cell(70,1);
       imagePerCategory = 20;
       categoryDelim = '-';
       
       lastClass = '';
       classItr = 1;
       for fileItr = 1:numel(fileList)
           [pathToImages, fileName, ext] = fileparts(fileList{fileItr});
           fileName
           imageClassIdx = strfind(fileName, categoryDelim);
           imageClass = fileName(1,1:(imageClassIdx(1)-1));
           
           if ~strcmp(imageClass, lastClass)
               lastClass = imageClass;
               classArr(classItr) = {imageClass};
               classItr = classItr + 1;
           end
       end
       for classItr = 1:70
           %% Assign some random images from each category as the vocabulary learning dataset.
           imageIdx = randsample(imagePerCategory, numberOfVocabImagesPerCategory);
           for fileItr = 1:numberOfVocabImagesPerCategory
               if exist([datasetFolder '/' classArr{classItr} '-0' num2str(imageIdx(fileItr)) '.gif'], 'file')
                    newImagePath = [vocabImageFolder '/' classArr{classItr} '-0' num2str(imageIdx(fileItr)) '.gif'];
               else
                    imagePath = [datasetFolder '/' classArr{classItr} '-' num2str(imageIdx(fileItr)) '.gif'];
                    newImagePath = [vocabImageFolder '/' classArr{classItr} '-' num2str(imageIdx(fileItr)) '.gif'];
               end
               copyfile(imagePath, newImagePath);
           end
       end
    end
    
end