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
%> @param vocabLevel The redundant vocabulary level to search for.
%> @param previousGraphLevel The input graph's path.
%> @param options Program options.
%>
%> @retval graphLevel Realizations of each composition in the vocabulary.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.12.2013
function [newLevel] = collectInstances(vocabLevel, redundantVocabLevel, previousGraphLevel, options, levelItr)
    %% Discover new level's subs.
    [~, newLevel, ~] = discoverSubs(vocabLevel, redundantVocabLevel, previousGraphLevel, options, true, levelItr);
end
