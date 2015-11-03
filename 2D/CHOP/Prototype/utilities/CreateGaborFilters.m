%> Name: CreateGaborFilters
%>
%> Description: This function creates basic Gabor filters to 
%> detect low-level features in the input images. 
%>
%> @param gaborOptions Options for creating gabor filters.
%>
%> @retval filters n filters concatenated to a matrix of size d x d x n.
%> @retval filterImages RGB (uint8) images of the filters.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 25.10.2015
function [ filters, filterImages ] = CreateGaborFilters( gaborOptions )
    filterType = 'gabor';
    %% Gabor filters: A reference Gabor filter is rotated options.numberOfGaborFilters times to obtain steerable feature filters.
    filterFolder = [pwd '/filters/'];
    if ~exist([filterFolder filterType], 'dir')
       mkdir([filterFolder filterType]); 
    end
    
   %% Create reference filter.
   refGaborFilt = gabor_fn(gaborOptions.gabor.sigma, gaborOptions.gabor.theta, ...
       gaborOptions.gabor.lambda, ...
       gaborOptions.gabor.psi, gaborOptions.gabor.gamma);
   dimSize = size(refGaborFilt);

   %% Pad or crop the sides of the gabor filter to set it to a fixed size of gaborFilterSize x gaborFilterSize
   newFilt = zeros(gaborOptions.gaborFilterSize);
   if dimSize(1) >= gaborOptions.gaborFilterSize
      lowerPad = ceil((dimSize(1)-gaborOptions.gaborFilterSize)/2);
      upperPad = (dimSize(1) - lowerPad) - gaborOptions.gaborFilterSize;
      refGaborFilt = refGaborFilt((lowerPad+1):(end-upperPad), :);
   end
   if dimSize(2) >= gaborOptions.gaborFilterSize
      leftPad = ceil((dimSize(2)-gaborOptions.gaborFilterSize)/2);
      rightPad = (dimSize(2) - leftPad) - gaborOptions.gaborFilterSize;
      refGaborFilt = refGaborFilt(:, (leftPad+1):(end-rightPad));
   end
   dimSize = size(refGaborFilt);
   lowerPad = ceil((gaborOptions.gaborFilterSize - dimSize(1))/2);
   upperPad = (gaborOptions.gaborFilterSize - lowerPad) - dimSize(1);
   leftPad = ceil((gaborOptions.gaborFilterSize - dimSize(2))/2);
   rightPad = (gaborOptions.gaborFilterSize - leftPad) - dimSize(2);
   newFilt((lowerPad+1):(end-upperPad), ((leftPad+1):(end-rightPad))) = refGaborFilt;
   refGaborFilt = newFilt;

   %% Get other filters by rotating the reference filter.
   filters = cell(gaborOptions.numberOfGaborFilters,1);
   for filtItr = 0:(gaborOptions.numberOfGaborFilters-1)
     % If the filter does not exist, create it.
     if ~exist([pwd '/filters/filt' num2str(filtItr+1) '.mat'], 'file') 
          theta = (-180/gaborOptions.numberOfGaborFilters) * filtItr;
          gaborFilt = imrotate(refGaborFilt, theta, 'bilinear', 'crop');

          % Normalize the filter so there are no discrepancies between
          % them.
          curConst = sum(sum(gaborFilt));
          gaborFilt = gaborFilt - curConst/(numel(refGaborFilt));
          filters{filtItr+1} = gaborFilt;
          gaborFilt = gaborFilt - min(min(gaborFilt));

          % Save filters
          save([filterFolder filterType '/filt' num2str(filtItr+1) '.mat'], 'gaborFilt');
          imwrite(gaborFilt, [filterFolder filterType '/filt' num2str(filtItr+1) '.png']);
          imwrite(gaborFilt>0, [filterFolder filterType '/filt' num2str(filtItr+1) 'Mask.png']);
     else
          load([filterFolder filterType '/filt' num2str(filtItr+1) '.mat']);
          filters{filtItr+1} = gaborFilt;
     end
   end
   filterImages = filters;
end

