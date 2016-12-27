%> Name: createFilters
%>
%> Description: This function extract basic Gabor/auto filter responses to 
%> detect low-level features in the input images. 
%>
%> @param options Program options.
%>
%> @retval filters n filters concatenated to a matrix of size d x d x n.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 03.12.2013 Response inhibition added.
%> Ver 1.2 on 12.01.2014 Comment changes to create unified code look.
function [ filters, filterImages ] = createFilters( options )
    %% Gabor filters: A reference Gabor filter is rotated options.numberOfGaborFilters times to obtain steerable feature filters.
    filterFolder = [options.currentFolder '/filters/'];
    if ~exist([filterFolder options.filterType], 'dir')
       mkdir([filterFolder options.filterType]); 
    end
    
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'realGabor')
        %% Create reference filter.
        refGaborFilt = gabor_fn(options.gabor.sigma, options.gabor.theta, ...
            options.gabor.lambda, ...
            options.gabor.psi, options.gabor.gamma);
        dimSize = size(refGaborFilt);

        %% Get other filters by rotating the reference filter.
        filters = cell(options.numberOfGaborFilters,1);
        for filtItr = 0:(options.numberOfGaborFilters-1)
          % If the filter does not exist, create it.
          if ~exist([options.currentFolder '/filters/filt' num2str(filtItr+1) '.mat'], 'file') 
               theta = (-180/options.numberOfGaborFilters) * filtItr;
               gaborFilt = imrotate(refGaborFilt, theta, 'bilinear');
               
               % CROP THE FILTER HERE.
               dimSize = size(gaborFilt);
               newFilt = zeros(options.gaborFilterSize);
               if dimSize(1) >= options.gaborFilterSize
                     lowerPad = ceil((dimSize(1)-options.gaborFilterSize)/2);
                     upperPad = (dimSize(1) - lowerPad) - options.gaborFilterSize;
                     gaborFilt = gaborFilt((lowerPad+1):(end-upperPad), :);
               end
               if dimSize(2) >= options.gaborFilterSize
                     leftPad = ceil((dimSize(2)-options.gaborFilterSize)/2);
                     rightPad = (dimSize(2) - leftPad) - options.gaborFilterSize;
                     gaborFilt = gaborFilt(:, (leftPad+1):(end-rightPad));
               end
               dimSize = size(gaborFilt);
               lowerPad = ceil((options.gaborFilterSize - dimSize(1))/2);
               upperPad = (options.gaborFilterSize - lowerPad) - dimSize(1);
               leftPad = ceil((options.gaborFilterSize - dimSize(2))/2);
               rightPad = (options.gaborFilterSize - leftPad) - dimSize(2);
               newFilt((lowerPad+1):(end-upperPad), ((leftPad+1):(end-rightPad))) = gaborFilt;
               gaborFilt = newFilt;

               % Normalize the filter so there are no discrepancies between
               % them.
               curConst = sum(sum(gaborFilt));
               gaborFilt = gaborFilt - curConst/(numel(refGaborFilt));
               filters{filtItr+1} = gaborFilt;
               gaborFilt = gaborFilt - min(min(gaborFilt));

               % Save filters
               save([filterFolder options.filterType '/filt' num2str(filtItr+1) '.mat'], 'gaborFilt');
               imwrite(gaborFilt, [filterFolder options.filterType '/filt' num2str(filtItr+1) '.png']);
               imwrite(gaborFilt>0, [filterFolder options.filterType '/filt' num2str(filtItr+1) 'Mask.png']);
          else
               load([filterFolder options.filterType '/filt' num2str(filtItr+1) '.mat']);
               filters{filtItr+1} = gaborFilt;
          end
        end
        filterImages = filters;
    elseif strcmp(options.filterType, 'auto')
        filters = cell(options.autoFilterCount,1);
        filterImages = cell(options.autoFilterCount,1);
        for filtItr = 0:(options.autoFilterCount-1)
            fileName = [filterFolder options.filterType '/filt' num2str(filtItr+1) '.mat'];
            if exist(fileName, 'file')
               load(fileName, 'gaborFilt');
               filters{filtItr+1} = gaborFilt; %#ok<NODEF>
            end
            
            fileName = [filterFolder options.filterType '/filt' num2str(filtItr+1) '.png'];
            if exist(fileName, 'file')
               gaborFilt = imread(fileName);
               filterImages{filtItr+1} = gaborFilt; %#ok<NODEF>
            end
        end
    else
        display('Error: Filter type not implemented (in createFilters.m).');
        filters = [];
        filterImages = [];
    end
end

