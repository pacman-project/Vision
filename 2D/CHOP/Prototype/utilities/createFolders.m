function [ ] = createFolders( options )

    % Create folder structures required later.
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
end

