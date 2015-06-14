% this function reads meshes

function [V, F, N] = meshRead(filename)


    if strfind(filename, '.obj')
        [V,N,F] = simpleObjReader(filename);
    end

    % TODO add codes for different file formats

end

