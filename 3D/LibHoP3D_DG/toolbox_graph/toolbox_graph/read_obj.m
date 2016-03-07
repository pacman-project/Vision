function [vertex,faces,normal] = read_obj(filename, sizeEst) % ,normal, DF, PI, sizeEst)  

% read_obj - load a .obj file.
%
%   [vertex,face,normal] = read_obj(filename);
%
%   faces    : list of facesangle elements
%   vertex  : node vertexinatates
%   normal : normal vector list
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin == 1
  sizeEst = 3000;
end

fid = fopen(filename);
if fid<0
    error(['Cannot open ' filename '.']);
end

frewind(fid);
a = fscanf(fid,'%c',1);
if strcmp(a, 'P')
    % This is the montreal neurological institute (MNI) specific ASCII facesangular mesh data structure.
    % For FreeSurfer software, a slightly different data input coding is
    % needed. It will be provided upon request.
    fscanf(fid,'%f',5);
    n_points=fscanf(fid,'%i',1);
    vertex=fscanf(fid,'%f',[3,n_points]);
    normal=fscanf(fid,'%f',[3,n_points]);
    n_faces=fscanf(fid,'%i',1);
    fscanf(fid,'%i',5+n_faces);
    faces=fscanf(fid,'%i',[3,n_faces])'+1;
    fclose(fid);
    return;
end

frewind(fid);

vertex = zeros(3, sizeEst);
faces = zeros(3, sizeEst);
normal = [];
% normal = zeros(3, sizeEst);
% DF = zeros(9, sizeEst);
% PI = zeros(2,sizeEst);
  
curF = 1;
curV = 1;
% curN = 1;
% curDF = 1;
% curPI = 1;
inds = [1,3,5];
temp = zeros(1, 6);

while 1
    s = fgetl(fid);
    if ~ischar(s), 
        break;
    end
    if ~isempty(s) && length(s)>2
        if strcmp(s(1:2), 'v ')
            % vertex
            vertex(:,curV) = sscanf(s(3:end), '%f %f %f');
            curV = curV + 1;
        end
        if strcmp(s(1:2), 'f ')
            % face
            
            % there are two possible notations
            if strfind(s, '//')
                 temp = sscanf(s(3:end), '%d//%d %d//%d %d//%d');
                 faces(:,curF) = temp(inds);
                 curF = curF + 1;
            else
                faces(:,curF) = sscanf(s(3:end), '%d %d %d');
                curF = curF + 1;
            end
        end
                %         if strcmp(s(1:2), 'vn')
                %             % normal
                %             normal(:,curN) = sscanf(s(3:end), '%f %f %f');
                %             curN = curN + 1;
                %         end
                %         if strcmp(s(1:2), 'pi')
                %             % part ID
                %             PI(:,curPI) = sscanf(s(3:end), '%d %d');
                %             curPI = curPI + 1;
                %         end
                %         if strcmp(s(1:2), 'df')
                %             % darboux frame
                %             DF(:,curDF) = sscanf(s(3:end), '%f %f %f %f %f %f %f %f %f');
                %             curDF = curDF + 1;
                %         end
    end
end

vertex = vertex(:, 1:curV - 1);
faces = faces(:, 1:curF - 1);
% PI = PI(:, 1:curPI - 1);
% normal = normal(:, 1:curN - 1);
% DF = DF(:, 1:curDF - 1);

fclose(fid);


