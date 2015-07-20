% this function reads meshes

function [V, F, N] = meshRead(filename)


    if strfind(filename, '.obj')
        [V, F, N] = read_obj(filename);
    end

    % TODO add codes for different file formats

end

