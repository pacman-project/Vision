%> Name: collectInstances
%>
%> Description: Given the vocabulary, this function detects realizations of
%> each composition in the vocabulary in the input graph. Although SUBDUE
%> hands out the realizations by default, sometimes they need to be collected
%> manually by using each separate definition (subgraph) as a pre-defined
%> substructure. SUBDUE essentially approximates a good solution, but a
%> comprehensive list of all instances are needed so that object graph is
%> formed correctly.
%>
%> @param vocabLevel The vocabulary level to search for.
%> @param graphFileName The input graph's path.
%> @param options Program options.
%>
%> @retval graphLevel Realizations of each composition in the vocabulary.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.12.2013
function [newLevel] = collectInstances(vocabLevel, graphFileName, options)
    numberOfPSFiles = numel(vocabLevel);
    newLevel = [];
    
    % Create temporary folder to put pre-defined subs in.
    tempFolder = tempname;
    mkdir(tempFolder);
    resultFileName = [tempFolder '/temp.txt'];
    
    % Put vocabLevel into a cell arr to make it compatible with pre-defined
    % file generator.
    tempVocabulary = cell(2,1);
    tempVocabulary(2) = {vocabLevel};
    preparePreDefinedFiles( tempFolder, tempVocabulary );
    
    %% Find realizations of each composition.
    for psItr = 1:numberOfPSFiles
        preDefinedFile = [tempFolder '/ps1/ps' num2str(psItr) '.g']; 
        
        if options.debug
           [~, psFileName, ~] = fileparts(preDefinedFile);
   %        display(['Collecting instances of ' psFileName '.']); 
        end
        %% Discover new level's subs.
        discoverSubs(graphFileName, resultFileName, options, options.currentFolder, preDefinedFile);
        [~, psLevel] = parseResultFile(resultFileName, options);
        % Assign instances correct labels.
        for instanceItr = 1:numel(psLevel)
           psLevel(instanceItr).labelId = psItr;
        end

        % Combine new level with the instances of this sub.
        if numel(psLevel)>0
            newLevel = [newLevel, psLevel];
        end
    end
    rmdir(tempFolder, 's');
end
