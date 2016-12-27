%> Name: getCombinations
%>
%> Description: Find part combinations up to maxSize. 
%>
%> @param 
%> 
%> @retval 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.11.2016
function [allCliques, labeledCliques] = getCombinations(prevLevel, rfSize, minSize, maxSize, beam, nsubs)
     imageIds = prevLevel(:,end);
     prevLevelPartCount = max(prevLevel(:,1)) + 1;
     [uniqueImageIds, uniqueImageIdx] = unique(imageIds, 'R2012a');
     uniqueImageIdx = cat(1, uniqueImageIdx, numel(imageIds)+1);
     
     % Save prev cliques and adjacency arrs for fast access.
     prevCliquesArr = cell(max(uniqueImageIds),1);
     adjMatrixArr = cell(max(uniqueImageIds),1);
     singleCliques = cell(max(uniqueImageIds),1);
     % Go through every image and find cliques.
     parfor imageItr = 1:numel(uniqueImageIds)
          nodes = prevLevel(uniqueImageIdx(imageItr):(uniqueImageIdx(imageItr+1)-1),:);
          nodePositions = single(nodes(:,2:3));
          
          % Calculate pairwise distances.
          distances = squareform(pdist(nodePositions, 'chebychev'));
          distances(triu(distances)>0) = 0;
          adjMatrix = distances<= (rfSize-1) & distances > 0;
          adjMatrix(1:(size(nodes,1)+1):size(nodes,1)^2) = 0;
          
          % Find cliques of nodes. 
          nodeCount = size(nodes,1);
          prevCliques = int32((1:nodeCount)');
          dummyArr = cat(2, prevCliques, zeros(nodeCount, maxSize-1, 'int32'));
          if minSize == 1
               cliques = dummyArr;
          else
               cliques = int32([]);
          end
          
          % Save data.
          prevCliquesArr{imageItr} = prevCliques;
          adjMatrixArr{imageItr} = adjMatrix;
          
          % Save cliques.
          cliques(cliques > 0) = cliques(cliques > 0) + (uniqueImageIdx(imageItr) - 1);
          singleCliques{imageItr} = cliques;
     end
     
     % Save single cliques.
     if minSize == 1
          singleCliques = cat(1, singleCliques{:});
          allCliques = singleCliques;
     else
          allCliques = [];
     end
     clear singleCliques;
          
     % Beam search to find new cliques.
     for sizeItr = 2:maxSize
          imageCliques = cell(max(uniqueImageIds),1);
          % Go through every image and find cliques.
          parfor imageItr = 1:numel(uniqueImageIds)
               % Retrieve data 
               prevCliques = prevCliquesArr{imageItr};
               adjMatrix = adjMatrixArr{imageItr};
               nodeCount = size(adjMatrix,1);
               dummyLabels = (1:nodeCount)';

               % Find cliques of nodes. 
               % Allocate space for new cliques.
               newCliques = cell(size(prevCliques,1),1);

               % Find nodes adjacent to all nodes in a clique.
               for newCliqueItr = 1:size(prevCliques,1)
                    relevantRow = prevCliques(newCliqueItr, :);
                    newAdjNodes = int32(dummyLabels(all(adjMatrix(:, relevantRow), 2)));
                    newCliques{newCliqueItr} = cat(2, relevantRow(ones(size(newAdjNodes,1),1), :), newAdjNodes);
               end

               % Pad and save new cliques.
               newCliques = cat(1, newCliques{:});
               prevCliques = newCliques;
               prevCliquesArr{imageItr} = prevCliques;
               newCliques = cat(2, newCliques, zeros(size(newCliques,1), maxSize - sizeItr, 1, 'int32'));
               if isempty(prevCliques)
                    continue;
               end
               % Save cliques.
               newCliques(newCliques > 0) = newCliques(newCliques > 0) + (uniqueImageIdx(imageItr) - 1);
               imageCliques{imageItr} = newCliques;
          end
          
          %% "Beam part of the search. We eliminate some of the previous cliques if they are not promising.
          % Two-stage elimination. We find unique subs, and eliminate them
          % based on two thresholds, one for carrying to the next layer
          % search, and one for saving them in output.
          imageCliqueCounts = cellfun(@(x) size(x,1), imageCliques);
          imageCliques = cat(1, imageCliques{:});
          
          % Determine unique parts.
          labeledCliques = imageCliques;
          labeledCliques(imageCliques>0) = prevLevel(imageCliques(imageCliques>0),1);
          labeledCliques = sort(labeledCliques,2);
          dummyArr = double(prevLevelPartCount).^((size(labeledCliques,2):-1:1) - 1);
          dummyArr = dummyArr(ones(size(labeledCliques,1),1),:);
          [~, ~, uniqueIdx] = unique(sum(dummyArr.*double(labeledCliques), 2));
          
          % Count how many times each combination has been encountered.
          cliqueCounts = hist(uniqueIdx, 1:max(uniqueIdx));
          sortedCliqueCounts = sort(cliqueCounts, 'descend');
          
          % Eliminate based on beam counts.
          if numel(sortedCliqueCounts) > beam
               beamThr = sortedCliqueCounts(beam);
               partsToExtend = (cliqueCounts >= beamThr)';
               cliquesToExtend = partsToExtend(uniqueIdx);
               idxArr = mat2cell(cliquesToExtend, imageCliqueCounts, 1);
               prevCliquesArr = cellfun(@(x,y) x(y, :), prevCliquesArr, idxArr, 'UniformOutput', false);
          end
          
          % Eliminate based on nsubs.
          if numel(sortedCliqueCounts) > nsubs
               nsubsThr = sortedCliqueCounts(nsubs);
               partsToSave = (cliqueCounts >= nsubsThr)';
               cliquesToSave = partsToSave(uniqueIdx);
               imageCliques = imageCliques(cliquesToSave, :);
          end
          
          % Save output.
          allCliques = cat(1, allCliques, imageCliques);
     end
     
     %% Create candidate parts.
     labeledCliques = allCliques;
     labeledCliques(allCliques>0) = prevLevel(allCliques(allCliques>0),1);
     [labeledCliques, sortIdx] = sort(labeledCliques,2);
     sortIdx = sortIdx';
     sortIdx = sortIdx(:);
     addedVals = size(labeledCliques,2) * (floor(0:1/size(labeledCliques,2):(size(labeledCliques,1)-1/size(labeledCliques,2))))';
     sortIdx = sortIdx + addedVals;
     allCliques = allCliques';
     allCliques(:) = allCliques(sortIdx);
     allCliques = allCliques';
     
     % Sort rows according to the parts.
     dummyArr = double(prevLevelPartCount).^((size(labeledCliques,2):-1:1) - 1);
     dummyArr = dummyArr(ones(size(labeledCliques,1),1),:);
     [~, sortIdx] = sort(sum(dummyArr.*double(labeledCliques), 2));
     labeledCliques = labeledCliques(sortIdx, :);
     allCliques = allCliques(sortIdx,:);
end
