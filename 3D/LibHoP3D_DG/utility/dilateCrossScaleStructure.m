function crossScaleStructureOut = dilateCrossScaleStructure(crossScaleStructure, F, iter, lenF)

    fringAll = ComputeFringDeep(F, iter, lenF);
    fringAll = fringAll{iter};
    crossScaleStructureOut = crossScaleStructure;
    
    for i = 1:lenF
        if crossScaleStructure(i) == 0 % there is no inference at this point
            fringCur = fringAll{i};
            crossScaleStructureOut(fringCur) = 0;
        end
    end
end