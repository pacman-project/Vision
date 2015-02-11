function [statisticsLayer, statisticsLayerSieved, statisticsLayerAggregated, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, ...
    partsLayer, fileForVisualizationLayer, partsLayerLoc, partsLayerEnt, partsSpecialSelected, partsLayerAll] =  getStandardFilePaths(root, dsN, nCl, aL)

    
    % files for the raw statistics
statisticsLayer{1} = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
statisticsLayer{3} = [root, 'statistics/statistics_3_', dsN, '_', nCl, '.mat'];
statisticsLayer{4} = [root, 'statistics/statistics_4_', dsN, '_', nCl, '.mat'];
statisticsLayer{5} = [root, 'statistics/statistics_5_', dsN, '_', nCl, '.mat'];
statisticsLayer{6} = [root, 'statistics/statistics_6_', dsN, '_', nCl, '.mat'];
statisticsLayer{7} = [root, 'statistics/statistics_7_', dsN, '_', nCl, '.mat'];
statisticsLayer{8} = [root, 'statistics/statistics_8_', dsN, '_', nCl, '.mat'];

% files for the sieved statistics
statisticsLayerSieved_Weak{1} = [root, 'statistics/statisticsSieved_Weak_1_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved_Weak{3} = [root, 'statistics/statisticsSieved_Weak_3_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved_Weak{4} = [root, 'statistics/statisticsSieved_Weak_4_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved_Weak{5} = [root, 'statistics/statisticsSieved_Weak_5_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved_Weak{6} = [root, 'statistics/statisticsSieved_Weak_6_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved_Weak{7} = [root, 'statistics/statisticsSieved_Weak_7_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved_Weak{8} = [root, 'statistics/statisticsSieved_Weak_8_', dsN, '_', nCl, '.mat'];

% files for the Aggregated statistics
statisticsLayerAggregated_Weak{1} = [root, 'statistics/statisticsAggregated_Weak_1_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated_Weak{3} = [root, 'statistics/statisticsAggregated_Weak_3_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated_Weak{4} = [root, 'statistics/statisticsAggregated_Weak_4_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated_Weak{5} = [root, 'statistics/statisticsAggregated_Weak_5_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated_Weak{6} = [root, 'statistics/statisticsAggregated_Weak_6_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated_Weak{7} = [root, 'statistics/statisticsAggregated_Weak_7_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated_Weak{8} = [root, 'statistics/statisticsAggregated_Weak_8_', dsN, '_', nCl, '.mat'];

% files for the sieved statistics
statisticsLayerSieved{1} = [root, 'statistics/statisticsSieved_1_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved{3} = [root, 'statistics/statisticsSieved_3_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved{4} = [root, 'statistics/statisticsSieved_4_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved{5} = [root, 'statistics/statisticsSieved_5_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved{6} = [root, 'statistics/statisticsSieved_6_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved{7} = [root, 'statistics/statisticsSieved_7_', dsN, '_', nCl, '.mat'];
statisticsLayerSieved{8} = [root, 'statistics/statisticsSieved_8_', dsN, '_', nCl, '.mat'];

% files for the Aggregated statistics
statisticsLayerAggregated{1} = [root, 'statistics/statisticsAggregated_1_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated{3} = [root, 'statistics/statisticsAggregated_3_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated{4} = [root, 'statistics/statisticsAggregated_4_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated{5} = [root, 'statistics/statisticsAggregated_5_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated{6} = [root, 'statistics/statisticsAggregated_6_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated{7} = [root, 'statistics/statisticsAggregated_7_', dsN, '_', nCl, '.mat'];
statisticsLayerAggregated{8} = [root, 'statistics/statisticsAggregated_8_', dsN, '_', nCl, '.mat'];

% files for part selection results
partsLayer{3} = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayer{4} = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayer{5} = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayer{6} = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayer{7} = [root, 'statistics/partsSelectionResults_7_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayer{8} = [root, 'statistics/partsSelectionResults_8_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for part selection results
partsLayerAll{3} = [root, 'statistics/partsSelectionResultsAll_3_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerAll{4} = [root, 'statistics/partsSelectionResultsAll_4_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerAll{5} = [root, 'statistics/partsSelectionResultsAll_5_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerAll{6} = [root, 'statistics/partsSelectionResultsAll_6_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerAll{7} = [root, 'statistics/partsSelectionResultsAll_7_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerAll{8} = [root, 'statistics/partsSelectionResultsAll_8_', dsN, '_', nCl, '_a', aL, '.mat'];

partsLayerLoc{3} = [root, 'statistics/partsSelectionResultsLoc_3_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerLoc{4} = [root, 'statistics/partsSelectionResultsLoc_4_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerLoc{5} = [root, 'statistics/partsSelectionResultsLoc_5_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerLoc{6} = [root, 'statistics/partsSelectionResultsLoc_6_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerLoc{7} = [root, 'statistics/partsSelectionResultsLoc_7_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerLoc{8} = [root, 'statistics/partsSelectionResultsLoc_8_', dsN, '_', nCl, '_a', aL, '.mat'];

partsLayerEnt{3} = [root, 'statistics/partsSelectionResultsEnt_3_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerEnt{4} = [root, 'statistics/partsSelectionResultsEnt_4_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerEnt{5} = [root, 'statistics/partsSelectionResultsEnt_5_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerEnt{6} = [root, 'statistics/partsSelectionResultsEnt_6_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerEnt{7} = [root, 'statistics/partsSelectionResultsEnt_7_', dsN, '_', nCl, '_a', aL, '.mat'];
partsLayerEnt{8} = [root, 'statistics/partsSelectionResultsEnt_8_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for layer visualization
fileForVisualizationLayer{2} = [];
fileForVisualizationLayer{3} = [root, 'statistics/fileForVisualization_3_', dsN, '_', nCl, '.mat'];
fileForVisualizationLayer{4} = [root, 'statistics/fileForVisualization_4_', dsN, '_', nCl, '.mat'];
fileForVisualizationLayer{5} = [root, 'statistics/fileForVisualization_5_', dsN, '_', nCl, '.mat'];
fileForVisualizationLayer{6} = [root, 'statistics/fileForVisualization_6_', dsN, '_', nCl, '.mat'];
fileForVisualizationLayer{7} = [root, 'statistics/fileForVisualization_7_', dsN, '_', nCl, '.mat'];
fileForVisualizationLayer{8} = [root, 'statistics/fileForVisualization_8_', dsN, '_', nCl, '.mat'];


% parts selected based on some special criteria
partsSpecialSelected{3} = [root, 'statistics/partsSelectionResultsSpecial_3_', dsN, '_', nCl, '_a', aL, '.mat'];
partsSpecialSelected{4} = [root, 'statistics/partsSelectionResultsSpecial_4_', dsN, '_', nCl, '_a', aL, '.mat'];
partsSpecialSelected{5} = [root, 'statistics/partsSelectionResultsSpecial_5_', dsN, '_', nCl, '_a', aL, '.mat'];
partsSpecialSelected{6} = [root, 'statistics/partsSelectionResultsSpecial_6_', dsN, '_', nCl, '_a', aL, '.mat'];
partsSpecialSelected{7} = [root, 'statistics/partsSelectionResultsSpecial_7_', dsN, '_', nCl, '_a', aL, '.mat'];
partsSpecialSelected{8} = [root, 'statistics/partsSelectionResultsSpecial_8_', dsN, '_', nCl, '_a', aL, '.mat'];

end

