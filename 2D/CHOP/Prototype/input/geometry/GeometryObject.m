%> Class Name: GeometryObject
%>
%> Description: Contains centerCoords and edges that together define a geometric
%> shape. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.10.2014
classdef GeometryObject
    properties
        angle;
        translation;
        scale;
        nodes;
        centerCoords;
        edges;
        objectType;
        imageSize;
    end
    methods
        % Constructor class.
        function obj = GeometryObject(objectType, imageSize)
           obj = obj.initialize(objectType, imageSize); 
        end
        
        % Initialization of geometric elements. Can be called again to
        % revert to factory settings.
        function obj = initialize(obj, objectType, imageSize)
            obj.angle = 0;
            obj.translation = [0,0];
            obj.scale = 1;
            obj.imageSize = [imageSize, imageSize];
            obj.centerCoords = round([imageSize/2, imageSize/2]);
            switch objectType
                case 'square'
                    obj.nodes = [-0.5, -0.5; ...
                             -0.5, 0.5; ...
                             0.5, -0.5; ...
                             0.5, 0.5];
                         obj.edges = [1 2; 1 3; 2 4; 3 4];
                case 'triangle'
                    obj.nodes = [0, -0.67; ...
                                 -0.57, 0.33; ...
                                 0.57, 0.33];
                    obj.edges = [1 2; 2 3; 3 1];
                case 'star'
                    obj.nodes = [-0.5, -0.25; ...
                                 0, -0.55; ...
                                 0.5, -0.25; ...
                                 0.33, 0.45; ...
                                 -0.33, 0.45];
                    obj.edges = [1 3; 1 4; 2 4; 2 5; 3 5]; 
            end
            obj.centerCoords = [0, 0];
            obj.objectType = objectType;
        end
        
        % The render function applies current transformations, and returns
        % the points which are to be rendered on the image.
        function img = render(obj)
           img = zeros(obj.imageSize, 'uint8');
           if strcmp(obj.objectType, 'circle')
               % If the object is a circle, all we need is to apply
               % translation & scaling, and then find the points in 2D
               % image to be marked as 1. 
               coords = obj.centerCoords;
               coords = coords + obj.translation;
               coords = coords + repmat(round(obj.imageSize/2), size(coords,1), 1);
               img = MidpointCircle(img, round(obj.scale/2), ...
                   coords(:,1), coords(:,2), 255);
           else
               coords = obj.nodes;
               
               % Apply rotation around origin.
               R = [cosd(obj.angle) -sind(obj.angle); sind(obj.angle) cosd(obj.angle)];
               coords = (R * coords')';
               
               % Scaling.
               coords = coords * obj.scale;
               
               % Translation.
               coords = coords + repmat(obj.translation, size(coords,1), 1);
               coords = coords + repmat(round(obj.imageSize/2), size(coords,1), 1);
               
               % Draw lines and fill them in.
               drawnLinesStart = coords(obj.edges(:,1),:);
               drawnLinesEnd = coords(obj.edges(:,2),:);
               assignedIdx = drawline(drawnLinesStart, drawnLinesEnd, obj.imageSize);
               img(assignedIdx) = 255;
           end
 %          img = imfill(img, 'holes');
        end
    end
end



















































