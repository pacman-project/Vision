function [] = runAnalysis(datasetName, experimentStr)

    % First, we run the experiment.
    Experiment(datasetName, '.png');

    % Run node analysis.
    runDiscriminativeAnalysis(datasetName);
    plotDiscriminativeAnalysis(datasetName);
    
    % Copy all relevant plots and data to experiments folder.
     expFolder = [pwd '/experiments/' datasetName '/' experimentStr];
    if ~exist(expFolder, 'dir')
       mkdir(expFolder); 
    end

    % First, copy plots to this folder. 
    copyfile([pwd '/categorization/analysis/' datasetName '/plots'], [expFolder '/plots']);
    
    % Get likelihood calculations.
    copyfile([pwd '/output/' datasetName '/logLikelihood.mat'], [expFolder '/logLikelihood.mat']); 
    
    % Next, get image reconstructions.
    copyfile([pwd '/output/' datasetName '/reconstruction'], [expFolder '/reconstruction']);
    
    % Finally, change the name of the three output folders in order to
    % preserve the results. 
    if exist([pwd '/output/' datasetName '_' experimentStr], 'dir')
       rmdir([pwd '/output/' datasetName '_' experimentStr], 's');
    end
    if exist([pwd '/output/' datasetName], 'dir')
        movefile([pwd '/output/' datasetName], [pwd '/output/' datasetName '_' experimentStr]);
    end
    if exist([pwd '/debug/' datasetName '_' experimentStr], 'dir')
       rmdir([pwd '/debug/' datasetName '_' experimentStr], 's');
    end
    if exist([pwd '/debug/' datasetName], 'dir')
        movefile([pwd '/debug/' datasetName], [pwd '/debug/' datasetName '_' experimentStr]);
    end
    copyfile([pwd '/parameters/SetParameters' datasetName '.m'], expFolder);
    copyfile([pwd '/SetParameters.m'], expFolder);
    if exist([pwd '/categorization/analysis/' datasetName '_' experimentStr], 'dir')
       rmdir([pwd '/categorization/analysis/' datasetName '_' experimentStr], 's'); 
    end
    if exist([pwd '/categorization/analysis/' datasetName], 'dir')
        movefile([pwd '/categorization/analysis/' datasetName], [pwd '/categorization/analysis/' datasetName '_' experimentStr]);
    end
end