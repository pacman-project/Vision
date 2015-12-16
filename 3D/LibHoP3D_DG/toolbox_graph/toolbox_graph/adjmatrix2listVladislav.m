%   adjmatrix2list - convert from matrix adjacency representation
%       to list adjacency. The adjacency can be a weighted matrix, 
%       and unlinked vertices should have weight either 'Inf' or '<=0'.

function [adj_list] = adjmatrix2listVladislav(r,c)

% n = size(A,1);
% 
% parfor i=1:n
%     adj_listT{i} = find(A(i, :) > 0);
% end

%     b = num2cell(A, 2);
%     adj_list = cellfun(@find, b, 'UniformOutput', false);
% 
%     [r,c] = find(A);
%     [r, ids] = sort(r, 'ascend');
%     c = c(ids);
%     r = r';
%     c = c';
    
    curCell = 1;
    prevI = 1;
    for i = 1:length(r)
        if r(i) ~= curCell
            adj_list{curCell} = c(prevI:i-1);
%             dist_list{curCell} = dist(prevI:i-1);
            prevI = i;
            curCell = r(i);
        end
    end
    adj_list{curCell}  = c(prevI:length(c));
%     dist_list{curCell} = dist(prevI:length(c));
end

