%> Name: getNumberOfFilters
%>
%> Description: Return the correct number of filters.
%>
%> @param options Program options
%>
%> @retval numberOfFilters
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.01.2014
function [numberOfFilters] = getNumberOfFilters(options)
    if strcmp(options.filterType, 'gabor')
        numberOfFilters = options.numberOfGaborFilters;
    elseif strcmp(options.filterType, 'auto')
        numberOfFilters = options.autoFilterCount;
    else
        display('Error: Filter type not implemented (in getNumberOfFilters.m).');
        numberOfFilters = 0;
    end
end