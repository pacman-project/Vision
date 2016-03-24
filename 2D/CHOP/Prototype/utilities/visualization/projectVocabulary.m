%> Name: projectVocabulary
%>
%> Description: Given a vocabulary, this function projects all nodes to
%> their layer 1 nodes using modal reconstruction and saves the results.
%>
%> @param vocabulary Program vocabulary.
%>   
%> @retval vocabulary The vocabulary with modal experts attached to nodes.
%>  
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 23.03.2016
function vocabulary = projectVocabulary( vocabulary)

     nodeCounts = cellfun(@(x) numel(x), vocabulary);
     maxLevel = find(nodeCounts, 1, 'last');
     prevLevel = [];
     
     for levelItr = 1:maxLevel
          vocabLevel = vocabulary{levelItr};
          for vocabNodeItr = 1:numel(vocabLevel)
               nodes = [vocabNodeItr, 0, 0, levelItr];
               experts = projectNode(nodes, vocabulary, 'modal', true);
               
               if levelItr > 1
                    % Combine experts with the imaginations from the previous
                    % layer.
                    allExperts = cell(size(experts,1),1);
                    for expertItr = 1:size(experts,1)
                         newChildren = prevLevel(experts(expertItr,1)).modalExperts;
                         newChildren(:,2:3) = repmat(experts(expertItr,2:3), size(newChildren,1),1) + newChildren(:,2:3);
                         allExperts{expertItr} = newChildren;
                    end
                    experts = cat(1, allExperts{:});
               end
               
               % Write them back.
               vocabLevel(vocabNodeItr).modalExperts = experts;
          end
          vocabulary{levelItr} = vocabLevel;
          prevLevel = vocabLevel;
     end
end