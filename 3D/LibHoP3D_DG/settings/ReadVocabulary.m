% vocabulary is read to three global variables

function ReadVocabulary(layerID)

    global lenOut;
    global partsOut;
    global pairsAll;
    
    if layerID >= 3
        dd = load('Temp/layer3/partSelection3.mat');
        lenOut{3} = dd.nNClusters{3};
        partsOut{3} = dd.partsOut;
        pairsAll{3} = dd.pairsAll;
    end
    if layerID >= 4
        dd = load('Temp/layer4/partSelection4.mat');
        lenOut{4} = dd.nNClusters{4};
        partsOut{4} = dd.partsOut;
        pairsAll{4} = dd.pairsAll;
    end
    if layerID >= 5
        dd = load('Temp/Layer5/partSelection5.mat');
        lenOut{5} = dd.nNClusters{5};
        partsOut{5} = dd.partsOut;
        pairsAll{5} = dd.pairsAll;
    end
    if layerID >= 6
        dd = load('Temp/Layer6/partSelection6.mat');
        lenOut{6} = dd.nNClusters{6};
        partsOut{6} = dd.partsOut;
        pairsAll{6} = dd.pairsAll;
    end
    if layerID >= 7
        dd = load('Temp/Layer7/partSelection7.mat');
        lenOut{7} = dd.nNClusters{7};
        partsOut{7} = dd.partsOut;
        pairsAll{7} = dd.pairsAll;
    end
    if layerID >= 8
        dd = load('Temp/Layer8/partSelection8.mat');
        lenOut{8} = dd.nNClusters{8};
        partsOut{8} = dd.partsOut;
        pairsAll{8} = dd.pairsAll;
    end
    if layerID >= 9
        dd = load('Temp/Layer6/partSelection9.mat');
        lenOut{9} = dd.nNClusters{9};
        partsOut{9} = dd.partsOut;
        pairsAll{9} = dd.pairsAll;
   end
        if layerID >= 10
        dd = load('Temp/Layer10/partSelection10.mat');
        lenOut{10} = dd.nNClusters{10};
        partsOut{10} = dd.partsOut;
        pairsAll{10} = dd.pairsAll;
    end
end

