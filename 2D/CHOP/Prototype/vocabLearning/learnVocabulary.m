%> Name: learnVocabulary
%>
%> Description: Given the graph of first-level responses and their
%> relations, this function extracts a hierarchical vocabulary by detecting
%> structural patterns inherent in the input graph. The pattern discovery
%> phase is hierarchic, and continues until a the graph is compressed into
%> a single node.
%>
%> @param allNodes The nodes of the level-1 input graph.
%> @param allEdges The edges of the level-1 input graph.
%> @param graphFileName The input graph.
%> @param resultFileName The file which will contain Subdue's output at
%> each level.
%> @param options Program options
%>
%> @retval vocabulary The hierarchic vocabulary learnt from the data.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.11.2013
function [ vocabulary ] = learnVocabulary( allNodes, allEdges, graphFileName, ...
                                                            resultFileName,...
                                                            options)
    % Allocate space for levels
    vocabulary = cell(options.maxLevels,1);
    nsubs = options.numberOfFilters;
    
    for levelItr = 1:options.maxLevels
        % Set maximum number of subs at 2nd layer so later layers use the
        % same number.
        if levelItr == 2
            nsubs = options.subdue.nsubs;
        end
    
        level(nsubs) = struct('id', [], 'position', [], 'children', [], 'parent', []);

        %% Run SUBDUE on the graph for the first time to go from level 1 to level 2.
        % Set the search options 
        subdueOptions = ' ';
        if options.subdue.diverse
           subdueOptions = [subdueOptions '-diverse '];
        end
        if options.subdue.valuebased
           subdueOptions = [subdueOptions '-valuebased '];
        end
        if options.subdue.overlap
           subdueOptions = [subdueOptions '-overlap '];
        end

        subdueOptions = [subdueOptions ' -nsubs ' num2str(options.subdue.nsubs) ...
                            ' -minsize ' num2str(options.subdue.minSize) ...
                            ' -maxsize ' num2str(options.subdue.maxSize) ...
                            ' -beam ' num2str(options.subdue.beam) ...
                            ' -fileoutput 3 ' ...
                            ' -out ' resultFileName ' ' graphFileName];

        % Form the command and run subdue with specified options
        command = ['./subdue' subdueOptions];
        system(command);

        % Create the new node list here
        newNodes = getNewNodes(resultFileName, options);
        
        % Get the edges 
        
    
        % Clear level structure
        clear level;
    
    end
end


%> Name: getNewNodes
%>
%> Description: This function parses the SUBDUE's result file and extracts
%> the new nodes to be used in the new graph (i.e. next level).
%>
%> @param resultFileName The SUBDUE's output file.
%> @param options Program options.
%>
%> @retval newNodes The node list including an identifier for each new
%> node.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 27.11.2013
function[newNodes] = getNewNodes(resultFileName, options)
    % Read result file into a string
    resultFileString = fileread(resultFileName);
    resultStringLength = numel(resultFileString);
    
    % Learn number of instances
    instanceIdx = strfind(resultFileString, options.subdue.instanceIndicator);
    numberOfNewNodes = numel(instanceIdx);
    
    % Learn number of types
    labelIdx = strfind(resultFileString, options.subdue.labelIndicator);
    labelIdx = labelIdx + numel(options.subdue.labelIndicator);
    numberOfNewLabels = numel(labelIdx);
    labelIds = cell(numberOfNewLabels,1);
    
    % Try to read all labels within the frames specified with
    % options.subdue.maxLabelLength.
    safeIdx = resultStringLength - options.maxLabelLength;
    for labelItr = 1:numberOfNewLabels
       labelStartIdx = labelIdx(labelItr);
       if labelIdx(labelItr) > safeIdx
          labelEndIdx = resultStringLength; 
       else
          labelEndIdx = labelStartIdx + options.maxLabelLength;
       end
       
       label = textscan(resultFileString(1,labelStartIdx:labelEndIdx), '%s', 1);
       labelIds(labelItr) = label{1};
    end
    
    %% Go over each instance and get ther label ids
    newNodes = cell(numberOfNewNodes,3);
    currentLabelIdx = 1;
    labelIdx = [labelIdx, resultStringLength];
    for instanceItr = 1:numberOfNewNodes
        % Find the label this instance belongs to
        while ~((labelIdx(currentLabelIdx) < instanceIdx(instanceItr)) && ...
                (instanceIdx(instanceItr) < labelIdx(currentLabelIdx+1)))
            currentLabelIdx = currentLabelIdx+1;
        end
        
        % Get the list of nodes from previous level which connect to this
        % instance.
        
        
        
        % Find the coordinate of this node by taking the center of gravity
        % of the nodes contributing to this node from previous layer. 
        % This step can be changed.
        
        
        newNodes{instanceItr,1} = labelIds{currentLabelIdx};
    end
end
