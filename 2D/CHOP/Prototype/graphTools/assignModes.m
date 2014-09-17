%> Name: assignModes
%>
%> Description: Given samples, assigns modes to each sample.
%> This function can be modified to introduce new mode calculation methods.
%>
%> @param nodes row-wise data samples of the form [x1, x2, .. , xd; y1,
%> ...].
%> @param minSamplesPerMode Minimum samples per mode, the function tries to
%> match this number for every cluster, by reducing number of clusters.
%> @param maximumModes Max. number of modes allowed.
%> @retval classes classes assigned to samples ( numberOfSamples x 1
%> matrix)
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.01.2014
function [ classes ] = assignModes( samples, minSamplesPerMode, maximumModes )
   classes = [];
   if isempty(samples)
       return;
   elseif size(samples,1) < minSamplesPerMode * maximumModes
        % Not enough samples, lower expected number of clusters.
        numberOfClusters = floor(size(samples,1)/minSamplesPerMode);
        if numberOfClusters <= 1
            classes = ones(size(samples,1),1);
        else
            classes = mec(samples, 'c', numberOfClusters, 'kmeans_i', 3);
        end
   else
       % Enough samples, process data.
       classes = mec(samples, 'c', maximumModes, 'kmeans_i', 3);
   end
end