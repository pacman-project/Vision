% creates a list of faces
% nf - number of faces to be extracted


function facesId = ExtractAdjacentFaces(faceCentral, nf, fring)

    % Create a list of faces
    facesId = faceCentral;

    while length(facesId) < nf
        if length(facesId) == 1
            facesId = [facesId, fring{facesId}];
        else
            facesIdNew = facesId;
            for j = 1:length(facesId)
                facesIdNew = union(facesIdNew, fring{facesId(j)});
            end
            facesId = facesIdNew;
        end
    end

end



%     startPoint = 1;
%     % EXTRACT POINTS AND NORMALS THAT BELONG TO THESE FACES
%     
%     if ~isempty(numPoints)
%         for j = 1:length(facesId)
%             endPoint = startPoint + numPoints(facesId(j)) - 1;
%             points(:,  startPoint:endPoint) = squeeze(VAll(facesId(j), :, 1:numPoints(facesId(j))));
%             Normals(:, startPoint:endPoint) = squeeze(NAll(facesId(j), :, 1:numPoints(facesId(j))));
%             startPoint = endPoint + 1;
%         end
%     elseif ~isempty(PI)
%         for jj = 1:nf
%             ids = find(PI(:,2) == facesId(jj));
%             endPoint = startPoint + length(ids) - 1;
% 
%             points(startPoint:endPoint, :) = VAll(ids, :);
%             normals(startPoint:endPoint, :) = NAll(ids, :);
%             startPoint = endPoint + 1;
%         end
%     end
%     
%     points = points(:, 1:endPoint);
%     Normals = Normals(:, 1:endPoint);

