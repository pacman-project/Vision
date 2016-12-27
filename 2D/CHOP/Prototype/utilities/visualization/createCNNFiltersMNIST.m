function [ output_args ] = createCNNFiltersMNIST(  )
     datasetName = 'MNIST';
     load([pwd '/output/' datasetName '/distributions.mat'], 'vocabularyDistributions');     
     load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary');
     rfSize = 5;
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
               if filterItr == 1
                    combinedPosProbs = combinedPosProbs(2:(end-1), 2:(end-1), :);
               else
                    combinedPosProbs = combinedPosProbs(1:(end-1), 1:(end-1), :);
               end
            
               %% Downsample combined probabilities.
               n = 2^(filterItr-1);
               k = rfSize;
               B = kron(speye(k), ones(1,n));
               downsampledPosProbs = zeros(rfSize, rfSize, numel(children));
               for childItr = 1:numel(children)
                    assignedVals = combinedPosProbs(:,:,childItr);
                    assignedVals = B*assignedVals*B';
                    downsampledPosProbs(:,:,childItr) = assignedVals;
               end
               
               % If RF is smaller, we move on and eliminate a row and a
               % column.
               if levelRFSizes(filterItr) < rfSize
                    allSums = sum(downsampledPosProbs,3);
                    columnSums = sum(allSums, 1);
                    elimColumn = rfSize;
                    if columnSums(1) < columnSums(end)
                         elimColumn = 1;
                    end
                    rowSums = sum(allSums, 2);
                    elimRow = rfSize;
                    if rowSums(1) < rowSums(end)
                         elimRow = 1;
                    end
                    downsampledPosProbs = downsampledPosProbs(setdiff(1:rfSize, elimRow), setdiff(1:rfSize, elimColumn), :);
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

