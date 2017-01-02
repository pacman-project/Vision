function [  ] = createCNNFiltersMNIST(  )
     datasetName = 'MNIST';
     load([pwd '/output/' datasetName '/distributions.mat'], 'vocabularyDistributions');     
     load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary');
     levelRFSizes = [5,5,4];
     filterItr = 1;
     f = 0.01;
     chopFilters = cell(numel(vocabularyDistributions) - 1,1);
     
     for vocabLevelItr = 2:numel(vocabularyDistributions)
          prevVocabLevel = vocabulary{filterItr};
          prevVocabLevelLabels = [prevVocabLevel.label];
          vocabLevel = vocabulary{vocabLevelItr};
          vocabLevelDistributions = vocabularyDistributions{vocabLevelItr};
          
          randCNNFilters = randn(levelRFSizes(filterItr), levelRFSizes(filterItr), numel(vocabularyDistributions{filterItr}), numel(vocabLevelDistributions), 'single');
          cnnFilters = zeros(levelRFSizes(filterItr), levelRFSizes(filterItr), numel(vocabularyDistributions{filterItr}), numel(vocabLevelDistributions), 'single');
          
          %% Go through the nodes and convert them to CNN filters.
          for vocabNodeItr = 1:numel(vocabLevelDistributions)
               tempFilters = zeros(levelRFSizes(filterItr), levelRFSizes(filterItr), numel(vocabularyDistributions{filterItr}), 'single');
               children = vocabLevel(vocabNodeItr).children;
               combinedPosProbs = vocabLevelDistributions(vocabNodeItr).childrenPosDistributionProbs;
               combinedPosProbs = cellfun(@(x) full(x), combinedPosProbs, 'UniformOutput', false);
               combinedPosProbs = cat(3, combinedPosProbs{:});
            
               %% Downsample combined probabilities.
               n = 2^(filterItr-1);
               k = levelRFSizes(filterItr);
               B = kron(speye(k), ones(1,n));
               downsampledPosProbs = zeros(levelRFSizes(filterItr), levelRFSizes(filterItr), numel(children));
               for childItr = 1:numel(children)
                    assignedVals = combinedPosProbs(:,:,childItr);
                    assignedVals = B*assignedVals*B';
                    downsampledPosProbs(:,:,childItr) = assignedVals;
               end
               
               %% Finally, assign values.
               for childItr = 1:numel(children)
                    assignedVals = downsampledPosProbs(:,:,childItr);
                    assignedVals = assignedVals * 5;
                    assignedVals = rot90(assignedVals,2);
                    realChildren = find(prevVocabLevelLabels == children(childItr));   
                    
                    for realChildItr = 1:numel(realChildren)
                         tempFilters(:,:,realChildren(realChildItr)) = max(tempFilters(:,:,realChildren(realChildItr)), assignedVals);
                    end
               end
               
               %% Assign the filter back.
               cnnFilters(:,:,:,vocabNodeItr) = tempFilters;
          end
          cnnFilters(cnnFilters == 0) = randCNNFilters(cnnFilters == 0);
          cnnFilters = f * cnnFilters;
          chopFilters{filterItr} = cnnFilters;
          filterItr = filterItr + 1;
     end
     save([pwd '/CNN/' datasetName '.mat'], 'chopFilters');
end

