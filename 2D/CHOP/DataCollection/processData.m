%> Process the data collected with the two cameras so that only the object
%> is contained in each image, not its surroundings.
function [ ] = processData( )
    angleStepper = 5;
    classStr = dir([pwd '/data']);
    
    for classItr = 1:numel(classStr)
       % Eliminate '.', '..', and '.DS_Store'
       if isempty(strfind(classStr(classItr).name, '.'))
           % Get Objects
           objectStr = dir([pwd '/data/' classStr(classItr).name]);
           
           for objectItr = 1:numel(objectStr)
               if isempty(strfind(objectStr(objectItr).name, '.'))
                   
                    % Set bounds for cropped images.
                    if strfind(objectStr(objectItr).date, '16-Apr-2014')
                        xLow = 160;
                        xHigh = 359;
                        yLow = 188;
                        yHigh = 387;
                    else
                        xLow = 160;
                        xHigh = 359;
                        yLow = 180;
                        yHigh = 379;
                    end
                    
                    % Create output folder
                    outputFolder = [pwd '/processedData/' classStr(classItr).name '/' objectStr(objectItr).name];
                    if exist(outputFolder, 'dir')
                        rmdir(outputFolder, 's');
                    end
                    mkdir(outputFolder);
                    
                    % Read each view, and crop it.
                    for viewItr = 0:angleStepper:(360-angleStepper)
                        img = imread([pwd '/data/' classStr(classItr).name '/' objectStr(objectItr).name '/left/' ...
                            objectStr(objectItr).name '_r' num2str(viewItr) '_l.png']);
                        img = img(xLow:xHigh, yLow:yHigh, :);
                        imwrite(img, [outputFolder '/' objectStr(objectItr).name '_r' num2str(viewItr) '_l.png']);
                    end
               end
           end
           
       end
    end
end

