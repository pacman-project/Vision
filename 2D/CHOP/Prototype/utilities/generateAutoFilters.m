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
function [ ] = generateAutoFilters( datasetName, fileType )
    %% Step 1: Get program options and initialize variables.
    options = SetParameters(datasetName, true);
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    fileNames = fuf([datasetFolder '*' fileType], 1, 'detail');
    
    if ~exist([options.currentFolder '/filters/auto'], 'dir')
        mkdir([options.currentFolder '/filters/auto']);
    else
        rmdir([options.currentFolder '/filters/auto'], 's');
        mkdir([options.currentFolder '/filters/auto']);
    end
    if ~exist([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'file')
        if ~exist([options.currentFolder '/filters/vis/' datasetName], 'dir')
            mkdir([options.currentFolder '/filters/vis/' datasetName]);
        end
        %% Set initial variables for collecting data.
        lowOffset = floor((options.autoFilterSize-1)/2);
        highOffset = floor(options.autoFilterSize/2);

        %% Step 2: Collect random samples from each image to create our data samples.
        display('Collecting samples...');
        samplesPerImg = ceil(options.autoFilterPatchCount / numel(fileNames));
        img = imread(fileNames{1});
        samples = zeros(samplesPerImg * numel(fileNames), (options.autoFilterSize^2) * size(img,3));

        sampleOffset = 0;
        for fileItr = 1:numel(fileNames)
            img = imread(fileNames{fileItr});
            
            %% Downsample the image if it will be processed that way.
            if max(size(img)) > options.maxImageDim
               img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
            end

            %% Get a pre-defined number of random centers from this image.
            validXInterval = (lowOffset+1):(size(img,1)-highOffset);
            validYInterval = (lowOffset+1):(size(img,2)-highOffset);
            centers = [datasample(validXInterval, samplesPerImg)', datasample(validYInterval, samplesPerImg)']; 
            invalidSampleCount = 0;
            for sampleItr = 1:samplesPerImg
                sample = img((centers(sampleItr,1)-lowOffset):(centers(sampleItr,1)+highOffset), ...
                    (centers(sampleItr,2)-lowOffset):(centers(sampleItr,2)+highOffset), :);
        %        size(sample)
                newSample = zeros(1, numel(sample));
                sampleStartItr = 1;
                pixelsPerBand = size(sample,1) * size(sample,2) - 1;
                for bandItr = 1:size(img,3)
                   sampleInd = sample(:,:,bandItr);
                   newSample(sampleStartItr:(sampleStartItr+pixelsPerBand)) = sampleInd;
                   sampleStartItr = sampleStartItr + pixelsPerBand + 1;
                end
                if size(samples,2) == numel(newSample)
                    samples(sampleItr+sampleOffset, :) = newSample;
                else
                    invalidSampleCount = invalidSampleCount + 1;
                end
            end
            sampleOffset = sampleOffset + samplesPerImg - invalidSampleCount;
        end
        
        %% Now, we'll apply ZCA whitening to all samples.
        display('Applying whitening...');
        [Xwh, mu, invMat, whMat] = whiten(samples); %#ok<ASGLU,NASGU>

        %% Cluster whitened samples to get cluster centers.
        display('Clustering features...');
        opts = statset('MaxIter', 600);
        [~, C] = kmeans(Xwh, options.autoFilterCount, 'Start', 'cluster', ...
            'EmptyAction', 'Singleton', 'Replicates', 1, 'Options', opts);

        %% Save cluster centers, along with other info.
        save([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'C', 'Xwh', 'mu', 'invMat', 'whMat');
    else
        return;
    end
    numberOfFilters = size(C,1);
    
    %% Visualize the centers as a grid image.
    display('Visualizing features...');
    visX = ceil(sqrt(numberOfFilters));
    visY = ceil(numberOfFilters/visX);
    imageSize = [visX * options.autoFilterSize + (visX-1), ...
        visY * options.autoFilterSize + (visY-1), size(img,3)];
    finalImage = zeros(imageSize);
    filterItr = 1;
    for xItr = 1:visX
        for yItr = 1:visY
            gaborFilt = reshape(C(((xItr-1)*visY)+yItr, :), [options.autoFilterSize, options.autoFilterSize, size(img,3)]);
            finalImage(((xItr-1)*(options.autoFilterSize+1)+1):(xItr*(options.autoFilterSize+1)-1), ...
                ((yItr-1)*(options.autoFilterSize+1)+1):(yItr*(options.autoFilterSize+1)-1), :) = gaborFilt;
            printedGaborFilt = (gaborFilt - min(min(min(gaborFilt)))) / (max(max(max(gaborFilt))) - min(min(min(gaborFilt))));
            imwrite(printedGaborFilt, [options.currentFolder '/filters/auto/filt' num2str(filterItr) '.png']);
            save([options.currentFolder '/filters/auto/filt' num2str(filterItr) '.mat'], 'gaborFilt');
            filterItr = filterItr+1;
        end
    end
    finalImage = (finalImage - min(min(min(finalImage)))) / (max(max(max(finalImage))) - min(min(min(finalImage))));
    imwrite(finalImage, [options.currentFolder '/filters/vis/' datasetName '/C.png']);
    close(gcf);
end

