%> Name: createSimilarityMatrix
%>
%> Description: Create initial distance matrix for CHOP. 
%> This function defines the similarity of level 1 parts, where 0 is the
%> closest (itself) and 1 is farthest.
%> TODO: Change name to createDistanceMatrix.
%>
%> @param options Program options.
%>
%> @retval simMat Distance matrix of level 1 parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.05.2014
%> Ver 1.1 on 12.08.2014 Adding "auto" type similarity matrix.
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
        % Oh god.. Shakespeare would be proud.
    elseif strcmp(options.filterType, 'auto')
        filters = options.filters;
        simMat = zeros(numel(filters));
        for filtItr = 1:(numel(filters)-1)
            filter1 = filters{filtItr};
            for filtItr2 = (filtItr+1):numel(filters)
                filter2 = filters{filtItr2};
                distance = sqrt(sum(sum(sum((filter1-filter2).^2))));
                simMat(filtItr, filtItr2) = distance;
                simMat(filtItr2, filtItr) = distance;
            end
        end
        simMat = simMat/max(max(simMat));
    else
        simMat = [];
    end
end

