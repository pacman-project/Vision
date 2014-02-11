%> Name: preparePreDefinedFiles
%>
%> Description: Given the vocabulary, this function prepares a number of
%> files including compositions at each vocabulary level, each defined as a
%> small graph consisting of nodes from the previous level.
%>
%> @param vocabulary The vocabulary learned from the data.
%>
%> @param preDefinedFolder The folder under which the created files are
%> put.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 14.12.2013
function [] = preparePreDefinedFiles( preDefinedFolder, vocabulary )
    %% For each vocabulary level l_i where i>0, create the substructures from their definitions in vocabulary.
    for levelItr = 1:(numel(vocabulary)-1)
        levelPSFolder = [preDefinedFolder, '/ps' num2str(levelItr)];
        % If it does not exist, create pre-defined substructure folder.
        if ~exist(levelPSFolder, 'dir');
            mkdir(levelPSFolder);
        end
        vocabLevel = vocabulary{levelItr+1};
        for compItr = 1:numel(vocabLevel)
            fd = fopen([preDefinedFolder, '/ps' num2str(levelItr) '/ps' num2str(compItr) '.g'], 'w');
            fprintf(fd, 'PS\n');
            composition = vocabLevel(compItr);
            printGraphToFile(fd, composition.children', composition.adjInfo, 1);
            fclose(fd);
        end
    end
end

