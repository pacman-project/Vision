function [ depthPath, elPath ] = getPathToData(dataSetNumber, commonRoot)

    if dataSetNumber == 1
        depthPathBasic = [commonRoot, 'Input Data/AimShape/20_views_1Scale'];     % D:\Input Data\AimShape\1T  
    elseif dataSetNumber == 2
        depthPathBasic = [commonRoot,'Input Data/Washington/Wash-rgbd-dataset_02_scale'];
    elseif dataSetNumber == 3
        depthPathBasic = [commonRoot, 'Input Data/VladislavSTD/Vladislav_STD/depth'];     
    end
    
    depthPath{1} = depthPathBasic;
    depthPath{2} = depthPathBasic;
    depthPath{3} = depthPathBasic;
    depthPath{4} = depthPathBasic;
    depthPath{5} = [depthPathBasic, '_D1'];
    depthPath{6} = [depthPathBasic, '_D1'];
    depthPath{7} = [depthPathBasic, '_D2'];
    depthPath{8} = [depthPathBasic, '_D2'];

    for i = 1:8
        elPath{i} = [depthPathBasic, '_layer', num2str(i)];
    end
    
end

