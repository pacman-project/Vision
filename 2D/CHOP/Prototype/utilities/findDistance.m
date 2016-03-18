%> Name: findDistance
%>
%> Description: This function finds the distance between two parts using a
%> number of different techniques. The core method is done using
%> 'modal', which calculates pairwise distance by comparing modal
%> reconstructions of the parts. 
%>
%> @param muArr1
%> @param muArr2
%> @param varArr1
%> @param varArr2
%> @param distType
%> @param searchHalfSize
%>
%> @retval distance
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.01.2016
function [ minDistance ] = findDistance( muArr1, muArr2, distType, searchHalfSize, searchMultiplier )
     % Program variables.
     imgSize = size(muArr1);
     lowVals = floor((imgSize - 2*searchHalfSize * searchMultiplier)/2);
     highVals = imgSize - (2*searchHalfSize*searchMultiplier + 1) - lowVals;
     centerPoint = lowVals + 1 + searchHalfSize * searchMultiplier;
     
     % Create random center points to test the hypothesis.
     [sampleCenterX, sampleCenterY] = meshgrid((centerPoint(1)-searchHalfSize*searchMultiplier):searchMultiplier:(centerPoint(1)+searchHalfSize*searchMultiplier), ...
          (centerPoint(2)-searchHalfSize * searchMultiplier):searchMultiplier:(centerPoint(2)+searchHalfSize * searchMultiplier));
     sampleCenterX = sampleCenterX(:);
     sampleCenterY = sampleCenterY(:);
     
     % For every pair of centers, try 
     if strcmp(distType, 'modal')
          minDistance = Inf;
          for sampleItr = 1:size(sampleCenterX,1)
               sampleMuArr1 = muArr1((sampleCenterX(sampleItr) - lowVals(1)):(sampleCenterX(sampleItr) + highVals(1)), ...
                    (sampleCenterY(sampleItr) - lowVals(2)):(sampleCenterY(sampleItr) + highVals(2)));
               sampleMuArr2 = muArr2((sampleCenterX(sampleItr) - lowVals(1)):(sampleCenterX(sampleItr) + highVals(1)), ...
                    (sampleCenterY(sampleItr) - lowVals(2)):(sampleCenterY(sampleItr) + highVals(2)));
               sampleVarArr1 = varArr1((sampleCenterX(sampleItr) - lowVals(1)):(sampleCenterX(sampleItr) + highVals(1)), ...
                    (sampleCenterY(sampleItr) - lowVals(2)):(sampleCenterY(sampleItr) + highVals(2)));
               sampleVarArr2 = varArr2((sampleCenterX(sampleItr) - lowVals(1)):(sampleCenterX(sampleItr) + highVals(1)), ...
                    (sampleCenterY(sampleItr) - lowVals(2)):(sampleCenterY(sampleItr) + highVals(2)));
             
               % Find total distance and return that.
               totalDistance = 0;
               for x = 1:size(sampleMuArr1,1)
                    for y = 1:size(sampleMuArr1,2)
                         dist1 = log(sampleVarArr2(x,y) / sampleVarArr1(x,y)) + (sampleVarArr1(x,y)^2 + (sampleMuArr1(x,y) - sampleMuArr2(x,y))^2) / (2 * sampleVarArr2(x,y)^2) - 1/2;       
                         dist2 = log(sampleVarArr1(x,y) / sampleVarArr2(x,y)) + (sampleVarArr2(x,y)^2 + (sampleMuArr2(x,y) - sampleMuArr1(x,y))^2) / (2 * sampleVarArr1(x,y)^2) - 1/2;        
                         dist = dist1 + dist2;
                         totalDistance = totalDistance + dist;
                    end
               end
               if totalDistance < minDistance
                    minDistance = totalDistance;
               end
          end
     end
end

