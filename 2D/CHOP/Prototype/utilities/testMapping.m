function [ testSuccess ] = testMapping( vocabLevel, graphLevel, newDistanceMatrix, previousGraphLevel )
     testSuccess = true;
     prevAdjInfo = {previousGraphLevel.adjInfo};
     prevLabelIds = [previousGraphLevel.labelId];
     wrongCounts = 0;
     matrixSize = size(newDistanceMatrix);
     
     faultIdx = zeros(numel(graphLevel),1) > 0;
     for graphLevelItr = 1:numel(graphLevel)
          children = graphLevel(graphLevelItr).children;
          orderedChildren = children;
          
          % Check for node labels.
          childrenTrueLabels = vocabLevel(graphLevel(graphLevelItr).labelId).children;
          ind = sub2ind(matrixSize, childrenTrueLabels, prevLabelIds(orderedChildren));
          if nnz(newDistanceMatrix(ind)) > 0
               faultIdx(graphLevelItr) = 1;
               wrongCounts = wrongCounts + 1;
          end
          
          % Check for edge labels.
          edges = prevAdjInfo{orderedChildren(1)};
          for edgeItr = 2:numel(children)
               edgeLabel = edges(edges(:,2) == orderedChildren(edgeItr), 3);
               if edgeLabel ~= vocabLevel(graphLevel(graphLevelItr).labelId).adjInfo(edgeItr -1, 3)
                    faultIdx(graphLevelItr) = 1;
                    wrongCounts = wrongCounts + 1;
               end
          end
     end
     
     % Finally, if the count of the faulty mappings are more than 1, return
     % false.
     if wrongCounts > 0
          testSuccess =false;
     end
end

