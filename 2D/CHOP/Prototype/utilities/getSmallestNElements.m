%> Name: getSmallestNElements
%>
%> Description: Given a vector, this function returns smallest n elements
%> of it.
%>
%> @param vect Data vector.
%> @param n Number of elements to be returned.
%> 
%> @retval idx Linear indices of elements returned.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.12.2013
function [idx] = getSmallestNElements(vect, n)
    [~, sortedIdx] = sort(vect);
    idx = sortedIdx(1:n);
    restVect = vect(sortedIdx((n+1):end));
    % Add nodes which qualify, since same-valued nodes are already in idx.
    % Number of such nodes cannot exceed n.
    restOfNodes = (find(restVect == vect(idx(end)))+n);
    if numel(restOfNodes)>n
       restOfNodes = restOfNodes(1:n); 
    end
    idx = [idx; restOfNodes];
    idx = unique(idx, 'stable');
end