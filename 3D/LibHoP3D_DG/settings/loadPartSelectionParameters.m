function [iterations, lenSelected, fieldSize, maxRelDepth, quant, partsCoverArea, numSimilar, threshPair, sieve_thresh] = loadPartSelectionParameters(dataSetNumber)

    iterations{3} = 150;
    iterations{4} = 300;
    iterations{5} = 400;
    iterations{6} = 500;
    iterations{7} = 600;
    iterations{8} = 700;
    
    lenSelected{3} = 14000;
    lenSelected{4} = 2000;
    lenSelected{5} = 500;
    lenSelected{6} = 600;
    lenSelected{7} = 800;
    lenSelected{8} = 900;
    
    fieldSize{3} = [17, 7,  10]; 
    fieldSize{4} = [17, 17, 10]; 
    fieldSize{5} = [53, 17, 53];
    fieldSize{6} = [53, 53, 53];
    fieldSize{7} = [161, 53, 501];
    fieldSize{8} = [161, 161, 701];
    
    maxRelDepth{3} = 200;
    maxRelDepth{4} = 300;
    maxRelDepth{5} = 2000;
    maxRelDepth{6} = 2000;
    maxRelDepth{7} = 3000;
    maxRelDepth{8} = 3000;
    
    quant{3} = 0.02;
    quant{4} = 0.02;
    quant{5} = 0.02;
    quant{6} = 0.02;
    quant{7} = 0.02;
    quant{8} = 0.02;
    
    partsCoverArea{2} = [5, 5,  71];
    partsCoverArea{3} = [15, 5,  71];
    partsCoverArea{4} = [15, 15, 101];
    partsCoverArea{5} = [15, 5,  71];
    partsCoverArea{6} = [15, 15, 101]; 
    partsCoverArea{7} = [17, 7,  71];
    partsCoverArea{8} = [17, 17, 101];
    
    numSimilar{3} = 7;
    numSimilar{4} = 7;
    numSimilar{5} = 13;
    numSimilar{6} = 19;
    numSimilar{7} = 23;
    numSimilar{8} = 25;
    
    sieveThreshMultiplier = 10^(-4);
    
    threshPair{3} = 4 * sieveThreshMultiplier; 
    sieve_thresh{3} = 2 * sieveThreshMultiplier;
    
    threshPair{4} = 3 * sieveThreshMultiplier; 
    sieve_thresh{4} = 2 * sieveThreshMultiplier;
    
    threshPair{5} = 1 * sieveThreshMultiplier; 
    sieve_thresh{5} = 0.5 * sieveThreshMultiplier;
    
    threshPair{6} = 0.2 * sieveThreshMultiplier; 
    sieve_thresh{6} = 0.1 * sieveThreshMultiplier;
    
    threshPair{7} = 0.01 * sieveThreshMultiplier; 
    sieve_thresh{7} = 0.003 * sieveThreshMultiplier;
    
    threshPair{8} = 0.002 * sieveThreshMultiplier; 
    sieve_thresh{8} = 0.001 * sieveThreshMultiplier;
    
    

end

