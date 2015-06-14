% this function returns a file with co-occurrence st

function [statisticsLayer, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, statisticsLayerSieved, statisticsLayerAggregated] = ...
                        GetStatisticsFiles(root, dsN, nCl, scAdder, layerID)
    
    % files for the sieved statistics
    lineLayer = num2str(layerID);
    
    statisticsLayer =                [root, 'statistics/statistics_',                lineLayer, '_', dsN, '_', nCl, '_', scAdder, '.mat'];
    statisticsLayerSieved_Weak =     [root, 'statistics/statisticsSieved_Weak_',     lineLayer, '_', dsN, '_', nCl, '_', scAdder, '.mat'];
    statisticsLayerAggregated_Weak = [root, 'statistics/statisticsAggregated_Weak_', lineLayer, '_', dsN, '_', nCl, '_', scAdder, '.mat'];
    statisticsLayerSieved =          [root, 'statistics/statisticsSieved_',          lineLayer, '_', dsN, '_', nCl, '_', scAdder, '.mat'];
    statisticsLayerAggregated =      [root, 'statistics/statisticsAggregated_',      lineLayer, '_', dsN, '_', nCl, '_', scAdder, '.mat'];
end

