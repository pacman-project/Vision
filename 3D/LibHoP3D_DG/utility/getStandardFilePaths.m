function [partsLayer, fileForVisualizationLayer, partsLayerLoc, partsLayerEnt, partsSpecialSelected, partsLayerAll, calibrationFile] =  getStandardFilePaths(root, dsN, nCl)


    % files for part selection results

    for i = 3:8
        partsLayer{i} =                  [root, 'statistics/partsSelectionResults_', num2str(i), '_', dsN, '_', nCl,  '.mat'];
        partsLayerAll{i} =            [root, 'statistics/partsSelectionResultsAll_', num2str(i), '_', dsN, '_', nCl,  '.mat'];
        partsLayerLoc{i} =            [root, 'statistics/partsSelectionResultsLoc_', num2str(i), '_', dsN, '_', nCl,  '.mat'];
        partsLayerEnt{i} =            [root, 'statistics/partsSelectionResultsEnt_', num2str(i), '_', dsN, '_', nCl,  '.mat'];
        fileForVisualizationLayer{i} =    [root, 'statistics/fileForVisualization_', num2str(i), '_', dsN, '_', nCl,  '.mat'];
        partsSpecialSelected{i} = [root, 'statistics/partsSelectionResultsSpecial_', num2str(i), '_', dsN, '_', nCl,  '.mat'];

    end

    % calibration file
    calibrationFile = [root, 'settings/calibration', dsN, '.mat'];

end


% % partsLayer{3} = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', '.mat'];
% % partsLayer{4} = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', '.mat'];
% % partsLayer{5} = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', '.mat'];
% % partsLayer{6} = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', '.mat'];
% % partsLayer{7} = [root, 'statistics/partsSelectionResults_7_', dsN, '_', nCl, '_a', '.mat'];
% % partsLayer{8} = [root, 'statistics/partsSelectionResults_8_', dsN, '_', nCl, '_a', '.mat'];
% 
% % files for part selection results
% partsLayerAll{3} = [root, 'statistics/partsSelectionResultsAll_3_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerAll{4} = [root, 'statistics/partsSelectionResultsAll_4_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerAll{5} = [root, 'statistics/partsSelectionResultsAll_5_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerAll{6} = [root, 'statistics/partsSelectionResultsAll_6_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerAll{7} = [root, 'statistics/partsSelectionResultsAll_7_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerAll{8} = [root, 'statistics/partsSelectionResultsAll_8_', dsN, '_', nCl, '_a', '.mat'];
% 
% partsLayerLoc{3} = [root, 'statistics/partsSelectionResultsLoc_3_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerLoc{4} = [root, 'statistics/partsSelectionResultsLoc_4_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerLoc{5} = [root, 'statistics/partsSelectionResultsLoc_5_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerLoc{6} = [root, 'statistics/partsSelectionResultsLoc_6_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerLoc{7} = [root, 'statistics/partsSelectionResultsLoc_7_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerLoc{8} = [root, 'statistics/partsSelectionResultsLoc_8_', dsN, '_', nCl, '_a', '.mat'];
% 
% partsLayerEnt{3} = [root, 'statistics/partsSelectionResultsEnt_3_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerEnt{4} = [root, 'statistics/partsSelectionResultsEnt_4_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerEnt{5} = [root, 'statistics/partsSelectionResultsEnt_5_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerEnt{6} = [root, 'statistics/partsSelectionResultsEnt_6_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerEnt{7} = [root, 'statistics/partsSelectionResultsEnt_7_', dsN, '_', nCl, '_a', '.mat'];
% partsLayerEnt{8} = [root, 'statistics/partsSelectionResultsEnt_8_', dsN, '_', nCl, '_a', '.mat'];
% 
% % files for layer visualization
% fileForVisualizationLayer{2} = [];
% fileForVisualizationLayer{3} = [root, 'statistics/fileForVisualization_3_', dsN, '_', nCl, '.mat'];
% fileForVisualizationLayer{4} = [root, 'statistics/fileForVisualization_4_', dsN, '_', nCl, '.mat'];
% fileForVisualizationLayer{5} = [root, 'statistics/fileForVisualization_5_', dsN, '_', nCl, '.mat'];
% fileForVisualizationLayer{6} = [root, 'statistics/fileForVisualization_6_', dsN, '_', nCl, '.mat'];
% fileForVisualizationLayer{7} = [root, 'statistics/fileForVisualization_7_', dsN, '_', nCl, '.mat'];
% fileForVisualizationLayer{8} = [root, 'statistics/fileForVisualization_8_', dsN, '_', nCl, '.mat'];
% 
% % parts selected based on some special criteria
% partsSpecialSelected{3} = [root, 'statistics/partsSelectionResultsSpecial_3_', dsN, '_', nCl, '_a', '.mat'];
% partsSpecialSelected{4} = [root, 'statistics/partsSelectionResultsSpecial_4_', dsN, '_', nCl, '_a', '.mat'];
% partsSpecialSelected{5} = [root, 'statistics/partsSelectionResultsSpecial_5_', dsN, '_', nCl, '_a', '.mat'];
% partsSpecialSelected{6} = [root, 'statistics/partsSelectionResultsSpecial_6_', dsN, '_', nCl, '_a', '.mat'];
% partsSpecialSelected{7} = [root, 'statistics/partsSelectionResultsSpecial_7_', dsN, '_', nCl, '_a', '.mat'];
% partsSpecialSelected{8} = [root, 'statistics/partsSelectionResultsSpecial_8_', dsN, '_', nCl, '_a', '.mat'];

