% produces list of elements given list of depths
% it aclually replaces depthPath with elPath for all files

function [list_els] = makeElList(list_depth, depthPath, elPath)
    
    lenF = length(list_depth);
    lenDPW = length(depthPath);
    
    parfor i = 1:lenF
        curStr = list_depth{i};
        fileName = curStr(lenDPW+1:end);
        list_els{i} = [elPath, fileName];  
    end

end

