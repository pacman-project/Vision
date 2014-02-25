%> Name: discoverSubs
%>
%> Description: This function discovers parts and their realizations with
%> graphs defined by vocabLevel and graphLevel, as well as their printed
%> versions in graphFileName. 
%>
%> @param vocabLevel If preDefinedSearch is 1, the compositions in this
%> vocabulary level are searched in graphLevel. If 0, simply ignored.
%> @param graphLevel The current object graphs' level.
%> @param options Program options.
%> @param currentFolder Path to the workspace folder.
%> @param preDefinedSearch If true, supervised search is run. 
%> If empty, unsupervised SUBDUE runs over graphLevel.
%> @param levelItr current level id.
%>
%> @retval vocabLevel If preDefinedSearch is 1, [] is returned. If it is 0,
%> best compositions discovered are returned.
%> @retval graphLevel The graph consisting of part realizations discovered.
%> It's just a dummy graph including realization labels and children.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.01.2014
%> 'self' type search added on 05.02.2014
function [vocabLevel, graphLevel] = discoverSubs( vocabLevel, graphLevel, oppositeModes, options, currentFolder, preDefinedSearch, levelItr)
    matlabpool close;
    startTime = tic;
    if ~preDefinedSearch
        display(['.... Discovering compositions in level ' num2str(levelItr) '.']); 
    end
    if strcmp(options.subdue.implementation, 'exe')
        %% Specify name of the graph files.
        graphFileName = [options.currentFolder '/graphs/' options.datasetName '/train/' options.datasetName '_' num2str(levelItr) '.g'];
        resultFileName = [options.outputFolder '/' options.datasetName '.txt'];
        
        %% Print the object graphs to a file.
        if strcmp(options.subdue.implementation, 'exe')
            printGraphLevel(graphFileName, graphLevel);
        end
        
        if preDefinedSearch

            numberOfPSFiles = numel(vocabLevel);
            newLevel = cell(numberOfPSFiles,1);

            % Create temporary folder to put pre-defined subs in.
            tempFolder = tempname;
            mkdir(tempFolder);

            % Put vocabLevel into a cell arr to make it compatible with pre-defined
            % file generator.
            tempVocabulary = cell(2,1);
            tempVocabulary(2) = {vocabLevel};
            preparePreDefinedFiles( tempFolder, tempVocabulary );

            %% Find realizations of each composition.
            for psItr = 1:numberOfPSFiles
                preDefinedFile = [tempFolder '/ps1/ps' num2str(psItr) '.g']; 

                %% Discover new level's subs.
                [~, psLevel] = runSubdueExec(graphFileName, resultFileName, options, currentFolder, preDefinedFile);
                % Assign instances correct labels.
                for instanceItr = 1:numel(psLevel)
                   psLevel(instanceItr).labelId = psItr;
                end

                % Combine new level with the instances of this sub.
                if numel(psLevel)>0
                    newLevel{psItr} =  psLevel;
                end
            end
            graphLevel = cat(2, newLevel{:});
            rmdir(tempFolder, 's');
        else
            [vocabLevel, graphLevel] = runSubdueExec(graphFileName, resultFileName, options, currentFolder, []);
        end
    elseif strcmp(options.subdue.implementation, 'self')
        if preDefinedSearch
            [vocabLevel, graphLevel] = runSubdue(vocabLevel, graphLevel, oppositeModes, options, true);
        else
            [vocabLevel, graphLevel] = runSubdue(vocabLevel, graphLevel, oppositeModes, options, false);
        end
    end
    
    display(['.... Time elapsed: ' num2str(toc(startTime)) ' secs.']);
    display(['.... ' options.subdue.implementation ' type discovery is used. Found ' ...
        num2str(numel(graphLevel)) ' instances of ' num2str(numel(vocabLevel)) ' compositions.']);
    matlabpool('open', options.numberOfThreads);
end

