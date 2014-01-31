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
    fileNames = fuf([datasetFolder '*.jpg'], 1, 'detail');
    img = imread(fileNames{1});
    
    if ~exist([options.currentFolder '/output/' datasetName '/C.mat'], 'file')
        %% Set initial variables for collecting data.
        lowOffset = floor((options.autoFilterSize-1)/2);
        highOffset = floor(options.autoFilterSize/2);

        %% Step 2: Collect random samples from each image to create our data samples.
        samplesPerImg = ceil(options.autoFilterPatchCount / numel(fileNames));
        img = imread(fileNames{1});
        samples = zeros(samplesPerImg * numel(fileNames), (options.autoFilterSize^2) * size(img,3));

        sampleOffset = 0;
        for fileItr = 1:numel(fileNames)
            img = imread(fileNames{fileItr});

            %% Get a pre-defined number of random centers from this image.
            validXInterval = (lowOffset+1):(size(img,1)-highOffset);
            validYInterval = (lowOffset+1):(size(img,2)-highOffset);
            centers = [datasample(validXInterval, samplesPerImg)', datasample(validYInterval, samplesPerImg)']; 
            for sampleItr = 1:samplesPerImg
                sample = img((centers(sampleItr,1)-lowOffset):(centers(sampleItr,1)+highOffset), ...
                    (centers(sampleItr,2)-lowOffset):(centers(sampleItr,2)+highOffset), :);
                samples(sampleItr+sampleOffset, :) = sample(:)';
            end
            sampleOffset = sampleOffset + samplesPerImg;
        end

        %% Now, we'll apply ZCA whitening to all samples.
        [Xwh, mu, invMat, whMat] = whiten(samples);


        %% Cluster whitened samples to get cluster centers.
        [~, C] = kmeans(Xwh, options.autoFilterCount, 'Start', 'cluster');
 %       [XwhLabeled] = mec(Xwh, 'c', options.autoFilterCount);

        %% Save cluster centers, along with other info.
        save([options.currentFolder '/output/' datasetName '/C.mat'], 'C', 'Xwh', 'mu', 'invMat', 'whMat');
    else
        load([options.currentFolder '/output/' datasetName '/C.mat']);
    end
    C = C*invMat + repmat(mu, [options.autoFilterCount,1]);
    C = uint8(C);
    
    %% Visualize the centers as a grid image.
    imageSize = [options.autoFilterVisX * options.autoFilterSize + (options.autoFilterVisX-1), ...
        options.autoFilterVisY * options.autoFilterSize + (options.autoFilterVisY-1), size(img,3)];
    finalImage = zeros(imageSize);
    for xItr = 1:options.autoFilterVisX
        for yItr = 1:options.autoFilterVisY
            finalImage(((xItr-1)*(options.autoFilterSize+1)+1):(xItr*(options.autoFilterSize+1)-1), ...
                ((yItr-1)*(options.autoFilterSize+1)+1):(yItr*(options.autoFilterSize+1)-1), :) = ...
                reshape(C(((xItr-1)*options.autoFilterVisY)+yItr, :), [options.autoFilterSize, options.autoFilterSize, size(img,3)]);
        end
    end
    imwrite(uint8(finalImage), [options.currentFolder '/output/' datasetName '/C.png']);
    
end

