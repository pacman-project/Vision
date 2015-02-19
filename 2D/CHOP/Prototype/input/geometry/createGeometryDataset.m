%> Name: createGeometryDataset
%>
%> Description: Creates 2D images of simple geometric shapes such as
%> square, triangle, circle, star after applying affine transformations in 2D
%> plane. The set of transformations considered includes translation,
%> scaling, and rotation. 
%>
%> @param none
%>
%> @retval none
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.02.2015
function [] = createGeometryDataset( )
    % Here, we define the objects. 
    rotationCount = 12;
    scaleCount = 4;
    scaleStep = 1.3;
    initialSize = 40;
    imageSize = 300;
    imageValidSize = 120;
    translateCount = 5;
    dilation = 1;
    noiseDev = 0.05;
    objects = {'square', 'triangle', 'star', 'circle'};
    
    rotationAngles = 0:ceil(360/rotationCount):360;
    imgItr = 1;
    for objectItr = 1:numel(objects)
        for angle = rotationAngles
            for scaleItr = 1:scaleCount
                % Find random translateCount suitable translations.
                availableCoords = round(rand(translateCount,2)*imageValidSize) - ...
                    round(imageValidSize)/2;
                for translateItr = 1:translateCount
                    object = GeometryObject(objects{objectItr}, imageSize);
                    
                    
                    
                    
                    % Calculate necessary transformations.
                    % Adding some gaussian noise to the transformation
                    % parameters so the objects do not exactly match.
                    currentScale = initialSize * scaleStep^(scaleItr-1);
                    object.angle = angle + 180 * noiseDev * randn(1,1);
                    object.scale = round(currentScale + currentScale * noiseDev * randn(1,1));
                    object.translation = availableCoords(translateItr,:);
                    
                    % Apply transformations to the object at hand, and then 
                    % render the corresponding image.
                    img = object.render();
                    img = imdilate(img, strel('disk', dilation));
                    
                    % Write the generated image to a file.
                    if ~exist(['./Geometry/vocab/' objects{objectItr}], 'dir')
                        mkdir(['./Geometry/vocab/' objects{objectItr}]);
                    end
                    imwrite(img, ['./Geometry/vocab/' objects{objectItr} '/' num2str(imgItr) '.png']);
                    
                    % Increment image iterator.
                    imgItr = imgItr + 1;
                end
            end
        end
    end
end
