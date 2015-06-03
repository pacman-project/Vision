strRootE = 'D:\3D\Input Data\WashingtonLayers0.1\2Layer';

is_subset = false;
subsetPercent = 100;
[list_el, model_id, category_id, len] = extractFileListWashingtonForClassification(strRootE, is_subset, subsetPercent);

prevCat = 0;
prevMod = 0;
for i = 1:len
    
    I = imread(list_el{i});
    
    curMod = model_id(i);
    curCat = category_id(i);
    
    if prevCat ~= curCat
        str = ['curCat = ', num2str(curCat)];
        disp(str); 
        prevCat = curCat;
    end
    
    if prevMod ~= curMod
        str = ['curMod = ', num2str(curMod)];
        disp(str); 
        prevMod = curMod;
    end
       

end