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
     
     for levelItr = 1:maxLevel
          vocabLevel = vocabulary{levelItr};
          
          if isempty(vocabLevel) 
               break;
          end
          
          if ~isempty(vocabLevel(1).modalExperts)
               continue;
          end
          for vocabNodeItr = 1:numel(vocabLevel)
               nodes = [vocabNodeItr, 0, 0, levelItr];
               experts = projectNode(nodes, vocabulary, 'modal');
               
               % Write them back.
               vocabLevel(vocabNodeItr).modalExperts = experts;
          end
          vocabulary{levelItr} = vocabLevel;
     end
end