% this function finds the id of the face given 3 ids of vertices

function [fID] = findFace(F, V1, V2, V3)
    
    ids1 = F(1, :) == V1 | F(1, :) == V2 | F(1, :) == V3;
    ids2 = F(2, :) == V1 | F(2, :) == V2 | F(2, :) == V3;
    ids3 = F(3, :) == V1 | F(3, :) == V2 | F(3, :) == V3;
    
    fID = find(ids1 & ids2  & ids3 );
end

