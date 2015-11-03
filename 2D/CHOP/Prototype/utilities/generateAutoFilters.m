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
    if ~strcmp(options.filterType, 'auto')
        return;
    end
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    
    if ~exist([options.currentFolder '/filters/auto'], 'dir')
        mkdir([options.currentFolder '/filters/auto']);
    else
        rmdir([options.currentFolder '/filters/auto'], 's');
        mkdir([options.currentFolder '/filters/auto']);
    end
    if ~exist([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'file')
        if strcmp(fileType, '.mat')
            load([options.currentFolder '/input/' datasetName '/Auto.mat']);
            numberOfImages = size(images,1);
            imgDim = size(images,4);
        else
            fileNames = fuf([datasetFolder '*' fileType], 1, 'detail');
            numberOfImages = numel(fileNames);
            img = imread(fileNames{1});
            imgDim = size(img,3);
        end
        if ~exist([options.currentFolder '/filters/vis/' datasetName], 'dir')
            mkdir([options.currentFolder '/filters/vis/' datasetName]);
        end
        %% Set initial variables for collecting data.
        lowOffset = floor((options.autoFilterSize-1)/2);
        highOffset = floor(options.autoFilterSize/2);

        %% Step 2: Collect random samples from each image to create our data samples.
        display('Collecting samples...');
        samplesPerImg = ceil(options.autoFilterPatchCount / numberOfImages);
        samples = zeros(samplesPerImg * numberOfImages, (options.autoFilterSize^2) * imgDim);

        sampleOffset = 0;
        for fileItr = 1:numberOfImages
            if strcmp(fileType, '.mat')
                img = squeeze(images(fileItr,:,:,:));
            else
                img = imread(fileNames{fileItr});
            end
            
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
        opts = statset('MaxIter', 200);
        [~, C] = kmeans(Xwh, options.autoFilterCount, 'Start', 'cluster', ...
            'EmptyAction', 'Singleton', 'Display', 'iter', 'Options', opts);

        %% Save cluster centers, along with other info.
        save([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'C', 'mu', 'invMat', 'whMat');
%        save([options.currentFolder '/filters/vis/' datasetName '/Xwh.mat'], 'Xwh');
    else
        load([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'C', 'mu', 'invMat', 'whMat');
    end
    numberOfFilters = size(C,1);
 
    %% Visualize the centers as a grid image.
    display('Visualizing features...');
    visX = ceil(sqrt(numberOfFilters));
    visY = ceil(numberOfFilters/visX);
    imageSize = [visX * options.autoFilterSize + (visX-1), ...
        visY * options.autoFilterSize + (visY-1), round(size(C,2)/(options.autoFilterSize^2))];
    for imgItr = 1:2
         finalImage = zeros(imageSize);
         filterItr = 1;
         if imgItr == 1
              visC = C *invMat + repmat(mu, size( C,1),1);
         else
              visC = C;
              visC = (visC - min(min(visC))) / (max(max(visC)) - min(min(visC)));
         end
         if max(max(visC)) > 255
              visC = visC / 256;
         end
         for xItr = 1:visX
             for yItr = 1:visY
                 if ((xItr-1)*visY)+yItr > numberOfFilters
                     continue;
                 end
                 gaborFilt = reshape(C(((xItr-1)*visY)+yItr, :), [options.autoFilterSize, options.autoFilterSize, imageSize(3)]); %#ok<NASGU>
                 printedGaborFilt = reshape(visC(((xItr-1)*visY)+yItr, :), [options.autoFilterSize, options.autoFilterSize, imageSize(3)]);     
                 if imgItr == 2
                      printedGaborFilt = (printedGaborFilt - min(min(min(printedGaborFilt)))) / (max(max(max(printedGaborFilt))) - min(min(min(printedGaborFilt))));
                      save([options.currentFolder '/filters/auto/filt' num2str(filterItr) '.mat'], 'gaborFilt');
                 else
                      imwrite(uint8(round(printedGaborFilt)), [options.currentFolder '/filters/auto/filt' num2str(filterItr) '.png']);
                 end
                 finalImage(((xItr-1)*(options.autoFilterSize+1)+1):(xItr*(options.autoFilterSize+1)-1), ...
                      ((yItr-1)*(options.autoFilterSize+1)+1):(yItr*(options.autoFilterSize+1)-1), :) = printedGaborFilt;
                 filterItr = filterItr+1;
             end
         end
         if imgItr == 1
               finalImage = uint8(round(finalImage));
               imwrite(finalImage, [options.currentFolder '/filters/vis/' datasetName '/C_Im.png']);
         else
   %            finalImage = (finalImage - min(min(min(finalImage)))) / (max(max(max(finalImage))) - min(min(min(finalImage))));
               imwrite(finalImage, [options.currentFolder '/filters/vis/' datasetName '/C_True.png']);
         end
    end
    
    close(gcf);
end

