function [] = runAnalysis(datasetName, experimentStr)

    % First, we run the experiment.
    Experiment(datasetName, '.png');
    
    % Copy all relevant plots and data to experiments folder.
     expFolder = [pwd '/experiments/' datasetName '/' experimentStr];
    if exist(expFolder, 'dir')
       rmdir(expFolder, 's'); 
    end
     mkdir(expFolder);
    
    % First, copy plots to this folder. 
    try
          copyfile([pwd '/categorization/analysis/' datasetName '/plots'], [expFolder '/plots']);
    end
    
    try
          copyfile([pwd '/output/' datasetName], [expFolder '/output']);
    end
    
    try
         copyfile([pwd '/debug/' datasetName '/*.png'], [expFolder '/debug/']);
    end
    
    try
         copyfile([pwd '/parameters/SetParameters' datasetName '.m'], expFolder);
    catch
         copyfile([pwd '/parameters/SetParametersCommon.m'], expFolder);
    end
    
    try
          copyfile([pwd '/outputNodes/' datasetName '/nodes.mat'], [expFolder '/nodes.mat']);
    end
    
    try
         copyfile([pwd '/workspaces/' datasetName], [expFolder '/workspaces']);
    end
    
    try
         mkdir([expFolder '/models']);
         copyfile([pwd '/models/' datasetName '*.mat'], [expFolder '/models']);
    end
         
    try
         copyfile([pwd '/logs/' datasetName '.txt'], [expFolder '/' datasetName '.txt']);
    end
    
    copyfile([pwd '/SetParameters.m'], expFolder);
end