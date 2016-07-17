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
    isTraining = options.isTraining;
    try
    if isTraining 
       tempFolder = tempname;
       mkdir(tempFolder);
       if exist(options.outputFolder, 'dir')
           if ispc
               system(['move ' options.outputFolder tempFolder]);
           else
              rmdir(options.outputFolder, 's');
           end
       end
       if exist(options.workspaceFolder, 'dir')
           rmdir(options.workspaceFolder, 's');
       end
       if exist(options.debugFolder, 'dir')           
           if ispc
               system(['move ' options.debugFolder tempFolder]);
           else
              rmdir(options.debugFolder, 's');
           end
       end
       if exist(options.CNNFolder, 'dir')         
           rmdir(options.CNNFolder, 's');
       end
    end
    catch %#ok<CTCH>
        display('Output folders could not be removed. Images/data from output folder which is still open can cause such permission problems.');
    end
    
    if ~exist(options.processedFolder,'dir')
       mkdir(options.processedFolder); 
    end
    if ~exist(options.processedGTFolder,'dir')
       mkdir(options.processedGTFolder); 
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
    if ~exist(options.CNNFolder, 'dir')
         mkdir(options.CNNFolder);
    end
    if exist([options.currentFolder '/Workspace.mat'], 'file')
         delete([options.currentFolder '/Workspace.mat']);
    end
end

