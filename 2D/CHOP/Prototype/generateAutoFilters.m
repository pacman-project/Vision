%> Name: generateAutoFilters
%>
%> Description: This function clusters ZCA-whitened random samples from the
%> dataset, and clusters them to get a number of cluster centers as our basic
%> features for the hierarchy. The resulting features are used to detect
%> level 1 features (basic building blocks) of the vocabulary, while
%> compositions are defined as their combinations.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param imageExtension The extension of the files to work on. Examples
%> include '.jpg', '.png', '_crop.png'...
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 13.01.2014
function [ ] = generateAutoFilters( datasetName )
    %% Step 1: Get program options and initialize variables.
    options = SetParameters(datasetName);
    datasetFolder = [options.currentFolder '/output/' datasetName '/original/'];
    fileNames = fuf([datasetFolder '*.png'], 1, 'detail');
    
    %% Step 2: Collect random samples from each image to create our data samples.
    samplesPerImg = options.autoFilterPatchCount / numel(fileNames);
    samples = samplesPerImg * numel(fileNames);
    
    sampleOffset = 0;
    for fileItr = numel(fileNames)
        img = imread(fileNames{fileItr});
        
        centers = 
        for sampleItr = samplesPerImg
            
        end
        sampleOffset = sampleOffset + samplesPerImg;
    end
end

