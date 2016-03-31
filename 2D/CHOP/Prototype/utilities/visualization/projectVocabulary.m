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
function vocabularyDistributions = projectVocabulary( vocabularyDistributions)
     nodeCounts = cellfun(@(x) numel(x), vocabularyDistributions);
     maxLevel = find(nodeCounts, 1, 'last');
     
     for levelItr = 1:maxLevel
          vocabLevelDistributions = vocabularyDistributions{levelItr};
          
          if isempty(vocabLevelDistributions) 
               break;
          end
          
          if ~isempty(vocabLevelDistributions(1).modalExperts)
               continue;
          end
          for vocabNodeItr = 1:numel(vocabLevelDistributions)
               nodes = [vocabNodeItr, 0, 0, levelItr];
               experts = projectNode(nodes, vocabularyDistributions, 'modal');
               
               % Write them back.
               vocabLevelDistributions(vocabNodeItr).modalExperts = experts;
          end
          vocabularyDistributions{levelItr} = vocabLevelDistributions;
     end
end