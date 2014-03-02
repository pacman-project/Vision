%> Name: assignModes
%>
%> Description: Given samples, assigns modes to each sample.
%> This function can be modified to introduce new mode calculation methods.
%>
%> @param nodes row-wise data samples of the form [x1, x2, .. , xd; y1,
%> ...].
%> @param maximumModes Max. number of modes allowed.
%> @retval classes classes assigned to samples ( numberOfSamples x 1
%> matrix)
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.01.2014
function [ classes ] = assignModes( samples, maximumModes )
   classes = [];
   if isempty(samples)
       return;
   elseif size(samples,1) < maximumModes
       if size(samples,1) > 1
            classes = mec(samples, 'c', size(samples,1), 'kmeans_i', 5);
       else
            classes = 1;
       end
   else
       % Enough samples, process data.
       classes = mec(samples, 'c', maximumModes, 'kmeans_i', 5);
   end
end

