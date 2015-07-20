% darFrames should be of size n*9 or 9*n

function write_FirstLayer(filename, vertex, Norm, darFrames, partIDs)

% write_off - write a mesh to an OBJ file
%
%   write_obj(filename, vertex, face, options)
%
%   vertex must be of size [n,3]
%   face must be of size [p,3]
%
%   Copyright (c) 2015 Vladislav Kramarev

    %% check all the variables  
    
    vertex = CheckDimension(vertex, 3);
    Norm = CheckDimension(Norm, 3);
    partIDs = CheckDimension(partIDs, 2);
    darFrames = CheckDimension(darFrames, 9);
    
    % check if the folder exists
    folderName = returnFolder( filename );
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end

    %% open the file and write 
    fid = fopen(filename,'wt');
    if( fid==-1 )
        error('Can''t open the file.');
        return;
    end

    object_name = filename(1:end-4);

    fprintf(fid, '# write_obj (c) 2015 Vladislav Kramarev\n');

    object_name = 'curobj';

    fprintf(fid, ['g\n# object ' object_name ' to come\n']);

    % vertex position
    fprintf(fid, '# %d vertex\n', size(vertex,1));
    fprintf(fid, 'v %f %f %f\n', vertex');
    
    % write part IDs
    fprintf(fid, '\n');
    fprintf(fid, '# %d part Ids\n', size(partIDs,1));
    fprintf(fid, 'pi %d %d\n', partIDs');

    % Normals directions
    fprintf(fid, '\n');
    fprintf(fid, '# %d normals\n', size(Norm,1));
    fprintf(fid, 'vn %f %f %f\n', Norm');
    
    % Darboux frames
    fprintf(fid, '\n');
    fprintf(fid, '# %d Darboux Frames\n', size(darFrames,1));
    fprintf(fid, 'df %f %f %f %f %f %f %f %f %f\n', darFrames');

    fclose(fid);
end

function V = CheckDimension(V, correctSize)
    if size(V,2)~=correctSize
        V = V';
    end
    if size(V,2)~=correctSize
        error('An array does not have the correct format!');
    end
end



