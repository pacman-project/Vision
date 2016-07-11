function [graphLevel, vocabLevel, vocabularyDistributions] = updateDataStructures(graphLevel, vocabLevel, vocabularyDistributions, levelItr, options)
   % Update real label ids in main graph.
   graphRealLabels = [graphLevel.realLabelId];
   [remainingComps, ~, IC] = unique(graphRealLabels);
   IC = num2cell(int32(IC));
   [graphLevel.realLabelId] = deal(IC{:});
   
   % Remove unnecessary vocab node images. 
   images = cell(numel(remainingComps),1);
   for vocabNodeItr = 1:numel(remainingComps)
       images(vocabNodeItr) = {imread([options.debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(remainingComps(vocabNodeItr)) '.png'])};
   end
   delete([options.debugFolder '/level' num2str(levelItr) '/modalProjection/*.png']);
   
   % Update mu images.
   load([options.debugFolder '/level' num2str(levelItr) '/modalProjection/muImgs.mat']);
   smallMuImgs = smallMuImgs(remainingComps, :, :); %#ok<NODEF,NASGU>
   save([options.debugFolder '/level' num2str(levelItr) '/modalProjection/muImgs.mat'], 'smallMuImgs');
   
   for vocabNodeItr = 1:numel(remainingComps)
       imwrite(images{vocabNodeItr}, [options.debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '.png']);
   end
   
   % Eliminate unused compositions from vocabulary.
   vocabLevel = vocabLevel(remainingComps);
   vocabLevelDistributions = vocabularyDistributions{levelItr};
   vocabLevelDistributions = vocabLevelDistributions(remainingComps);
   vocabularyDistributions{levelItr} = vocabLevelDistributions;
   
   % Update vocabulary labels.
   vocabLevelLabels = [vocabLevel.label];
   [~, ~, IC] = unique(vocabLevelLabels);
   IC = num2cell(int32(IC));
   [vocabLevel.label] = deal(IC{:});
   
   % Update or node labels for the main graph.
   vocabLabels = [vocabLevel.label];
   graphRealLabels = [graphLevel.realLabelId];
   newGraphLabels = vocabLabels(graphRealLabels);
   newGraphLabels = num2cell(int32(newGraphLabels));
   [graphLevel.labelId] = deal(newGraphLabels{:});
   
   clear graphRealLabels IC remainingComps newGraphLabels vocabLabels;
end