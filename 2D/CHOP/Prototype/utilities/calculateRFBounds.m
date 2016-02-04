%> Name: calculateRFBounds
%>
%> Description: Given the positions of nodes in a higher layer, this function 
%> calculates the boundaries of the Receptive Field that can be represented
%> by each node. The output format is [minX, minY, maxX, maxY, ...].
%>
%> @param positions Node positions in level levelItr.
%> @param levelItr Level ID. 
%> @param RFSize Receptive field size.
%> @param poolDim Pooling size.
%>
%> @retval bounds RF Boundaries in terms of original pixels 
%> (not taking stride into account).
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.02.2016
function [ bounds ] = calculateRFBounds( positions, levelItr, options, isPrecise )
     % Get program parameters.
     RFSize = options.receptiveFieldSize;
     poolDim = options.poolDim;
     if isPrecise
          halfRFSize= 0;
     else
          halfRFSize = floor(RFSize/2);
     end
     if strcmp(options.filterType, 'gabor')
          stride = options.gabor.stride;
          inhibitionRadius = options.gabor.inhibitionRadius;
     else
          stride = options.auto.stride;
          inhibitionRadius = options.auto.inhibitionRadius;
     end

     % Calculate positions.
     positions = double(positions);
     if levelItr == 1
          bounds = [positions-halfRFSize, positions+halfRFSize];
     else
          downsampleFactor = poolDim ^ (levelItr-1);
          minStart = ((positions-1)-(2*halfRFSize));
          minStartPrecise = minStart * downsampleFactor + 1 + halfRFSize;
          maxEnd =  (positions+(2*halfRFSize));
          maxEndPrecise = maxEnd * downsampleFactor - (1 + halfRFSize);
          bounds = [minStartPrecise, maxEndPrecise];
     end
     
     % If precise boundaries are wanted, we process these guys differently.
     if isPrecise
          downsampleFactor = stride * (inhibitionRadius+1);
          bounds(:,1:2) = (bounds(:,1:2)-1)*downsampleFactor+1;
          bounds(:,3:4) = (bounds(:,3:4))*downsampleFactor;
     end
end

