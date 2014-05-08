%> Name: createSimilarityMatrix
%>
%> Description: Create initial similarity matrix for CHOP. 
%> This function defines the similarity of level 1 parts, where 0 is the
%> closest (itself) and 1 is farthest.
%> TODO: Make it adaptive to different number of filters and everything.
%>
%> @param options Program options.
%>
%> @retval simMat Similarity matrix of level 1 parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.05.2014
function [ simMat ] = createSimilarityMatrix( options )
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
        %% TODO: Change (CHANGE) the following similarity matrix definition with something more 'mathematical'. 
        % Oh god.
        simMat = [0 1/3 2/3 1 2/3 1/3; ...
                1/3 0 1/3 2/3 1 2/3; ...
                2/3 1/3 0 1/3 2/3 1; ... 
                1 2/3 1/3 0 1/3 2/3; ...
                2/3 1 2/3 1/3 0 1/3; ...
                1/3 2/3 1 2/3 1/3 0];
        % Oh god..
    else
        simMat = [];
    end
end

