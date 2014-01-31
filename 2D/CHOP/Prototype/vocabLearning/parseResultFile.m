%> Name: parseResultFile
%>
%> Description: This function parses the SUBDUE's result file and extracts
%> the new nodes to be used in the new graph (i.e. next level).
%>
%> @param resultFileName The SUBDUE's output file.
%> @param options Program options.
%> @param nsubs Number of subs allowed at this level of vocabulary.
%>
%>
%> @retval newNodes The node list including an identifier, position and a 
%> list of children for each new node.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.12.2013
function[vocabLevel, graphLevel] = parseResultFile(resultFileName, options)

    % Read result file into a string
    resultFileString = fileread(resultFileName);
    resultFileStringLength = numel(resultFileString);
    
    % Learn number of instances
    instanceIdx = strfind(resultFileString, options.subdue.instanceIndicator);
    numberOfNewNodes = numel(instanceIdx);
    
    % Learn number of label types
    labelIdx = strfind(resultFileString, options.subdue.labelIndicator);
    labelIdx = labelIdx + numel(options.subdue.labelIndicator);
    numberOfNewLabels = numel(labelIdx);
    labelIds = cell(numberOfNewLabels,1);
    
    if numberOfNewLabels == 0
       vocabLevel = []; 
       graphLevel = [];
       return;
    end
    
    % Allocate space for current vocabulary level.
    vocabLevel(numberOfNewLabels) = struct('label', [], 'children', [], 'parents', [], 'adjInfo', [], 'mdlScore', []);
    % Allocate space for current graph level.
    graphLevel(numberOfNewNodes) = struct('labelId', [], 'imageId', [], 'position', [], 'children', [], 'parents', [], 'adjInfo', [], 'childrenAdjInfo', []);
    
    % Try to read all labels within the frames specified with
    % options.subdue.maxLabelLength.
    safeIdx = resultFileStringLength - options.maxLabelLength;
    for labelItr = 1:numberOfNewLabels
       labelStartIdx = labelIdx(labelItr);
       if labelIdx(labelItr) > safeIdx
          labelEndIdx = resultFileStringLength; 
       else
          labelEndIdx = labelStartIdx + options.maxLabelLength;
       end
       
       label = textscan(resultFileString(1,labelStartIdx:labelEndIdx), '%s', 1);
       labelIds(labelItr) = label{1};
    end
    
    %% Go over each sub and extract info
    % Each label refers to a different sub in result file.
    currentInstanceOffset = 0;
    labelIdx = [labelIdx, resultFileStringLength];
    possibleSubDefnEndIdx = sort([labelIdx, instanceIdx, resultFileStringLength]);
    
    for labelItr = 1:numberOfNewLabels
        % From the definition of this sub, get node information for
        % vocabulary.
        subDefnStartIdx = labelIdx(labelItr);
        
        % TODO: Possible optimization here, decreased efficiency
        [subDefnEndIdx] = find(possibleSubDefnEndIdx>subDefnStartIdx, 1, 'first');
        subDefnEndIdx = possibleSubDefnEndIdx(subDefnEndIdx);
        
        % Read sub label and score here.        
        defnString = resultFileString(1, subDefnStartIdx:subDefnEndIdx);
        newLineIdx = strfind(defnString, sprintf('\n'));
        subLabel = sscanf(defnString(1:newLineIdx(1)), '%s%*s');
        mdlScore = sscanf(defnString((newLineIdx(1)+numel(options.subdue.scoreIndicator)): ...
                                                newLineIdx(2)), '%lf%*s');
        vocabLevel(labelItr).label = subLabel;
        vocabLevel(labelItr).mdlScore = mdlScore;
        
        % Find lines referring to nodes and edges in the definition.
        defnNodeIdx = strfind(defnString, options.subdue.nodeIndicator);
        defnEdgeIdx = strfind(defnString, options.subdue.edgeIndicator);
        defnDirEdgeIdx = strfind(defnString, options.subdue.directedEdgeIndicator);
        defnEdgeIdx = [defnEdgeIdx, defnDirEdgeIdx];
        
        %% Read node/edge information.
        childrenList = zeros(1, numel(defnNodeIdx));
        for nodeItr = 1:numel(defnNodeIdx)
           lineEnd = find(newLineIdx > defnNodeIdx(nodeItr), 1, 'first');
           nodeString = defnString(defnNodeIdx(nodeItr):newLineIdx(lineEnd));
           nodeLabel = textscan(nodeString, '%*s %*d %s');
           nodeLabel = nodeLabel{1};
           nodeLabel = nodeLabel{1};
           if strfind(nodeLabel, options.subdue.subPrefix)
              nodeLabel = str2double(nodeLabel((numel(options.subdue.subPrefix)+1):end));
           else
              nodeLabel = str2double(nodeLabel);
           end
           childrenList(nodeItr) = nodeLabel;
        end
        % Assign nodes to sub information
        vocabLevel(labelItr).children = childrenList;
        
        defnEdges = zeros(numel(defnEdgeIdx),4);
        for edgeItr = 1:numel(defnEdgeIdx)
           lineEnd = find(newLineIdx > defnEdgeIdx(edgeItr), 1, 'first');
           edgeString = defnString(defnEdgeIdx(edgeItr):newLineIdx(lineEnd));
           [vars] = textscan(edgeString, '%s %d %d %d');
           defnEdges(edgeItr,1:3) = cell2mat(vars(2:end));
           
           % Also save directedness information of the edge.               
           edgeChar = vars{1};
           if ~isempty(strfind(options.subdue.edgeIndicator, edgeChar{1}))
             defnEdges(edgeItr,4) = 0;
           else 
             defnEdges(edgeItr,4) = 1;
           end
        end
        vocabLevel(labelItr).adjInfo = defnEdges;
        
        %% Vocabulary-related info has been read. Now, parsing instances.
        % Find the label this instance belongs to
        subInstanceIdx = instanceIdx(instanceIdx>labelIdx(labelItr) & instanceIdx<labelIdx(labelItr+1));
        subInstanceIdx = [subInstanceIdx, labelIdx(labelItr+1)];
        for instanceItr = 1:(numel(subInstanceIdx)-1)
            instanceString =resultFileString(1,subInstanceIdx(instanceItr):subInstanceIdx(instanceItr+1));
            newLineIdx = strfind(instanceString, sprintf('\n'));
            nodeIdx = strfind(instanceString, options.subdue.nodeIndicator);
            
            % Determine current instance number and assign its label.
            currentInstance = instanceItr + currentInstanceOffset;
            graphLevel(currentInstance).labelId = labelItr;
            
            %% Read instance nodes.
            childrenList = zeros(1, numel(nodeIdx));
            for nodeItr = 1:numel(nodeIdx)
               lineEnd = find(newLineIdx > nodeIdx(nodeItr), 1, 'first');
               nodeString = instanceString(nodeIdx(nodeItr):newLineIdx(lineEnd));
               nodeId = textscan(nodeString, '%*s %d %*s');
               childrenList(nodeItr) = nodeId{1};
            end
            
            %% Read instance edges.
            edgeIdx = [strfind(instanceString, options.subdue.edgeIndicator), ...
                strfind(instanceString, options.subdue.directedEdgeIndicator)];
            instanceEdges = zeros(numel(edgeIdx),4);
            for edgeItr = 1:numel(edgeIdx)
               lineEnd = find(newLineIdx > edgeIdx(edgeItr), 1, 'first');
               edgeString = instanceString(edgeIdx(edgeItr):newLineIdx(lineEnd));
               [vars] = textscan(edgeString, '%s %d %d %d');
               instanceEdges(edgeItr,1:3) = cell2mat(vars(2:end));

               % Also save directedness information of the edge.
               edgeChar = vars{1};
               if ~isempty(strfind(options.subdue.edgeIndicator, edgeChar{1}))
                 instanceEdges(edgeItr,4) = 0;
               else 
                 instanceEdges(edgeItr,4) = 1;
               end
            end
            
            % Assign nodes to sub information
            graphLevel(currentInstance).children = childrenList;
            graphLevel(currentInstance).childrenAdjInfo = instanceEdges;
        end
        
        currentInstanceOffset = currentInstanceOffset + (numel(subInstanceIdx) - 1);
    end
end