%> Name: createFolders
%>
%> Description: Create the data structure required by CHOP to work.
%>
%> @param options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 16.06.2014
function [ ] = createFolders( options )
    if ~exist(options.processedFolder,'dir')
       mkdir(options.processedFolder); 
    end
    if ~exist(options.processedGTFolder,'dir')
       mkdir(options.processedGTFolder); 
    end
    if ~exist(options.preDefinedFolder,'dir')
       mkdir(options.preDefinedFolder); 
    end
    if ~exist(options.testGraphFolder, 'dir')
       mkdir(options.testGraphFolder);
    end
    if ~exist(options.trainGraphFolder, 'dir')
       mkdir(options.trainGraphFolder);
    end
    if ~exist(options.testOutputFolder, 'dir')
       mkdir(options.testOutputFolder); 
    end
    if ~exist(options.testInferenceFolder, 'dir')
       mkdir(options.testInferenceFolder); 
    end
    if ~exist(options.smoothedFolder,'dir')
       mkdir(options.smoothedFolder); 
    end
    if ~exist(options.debugFolder,'dir')
       mkdir(options.debugFolder); 
    end
end

