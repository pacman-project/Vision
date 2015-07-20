function [ depthPath ] = getPathToData(dataSetNumber, commonRoot)


    if dataSetNumber == 1
        depthPathBasic = [commonRoot, 'Input Data/AimShape/20_views_1Scale'];     % D:\Input Data\AimShape\1T  
    elseif dataSetNumber == 2
        depthPathBasic = [commonRoot,'Input data\Washington\Wash-rgbd-dataset_01_scale'];
    elseif dataSetNumber == 3
        depthPathBasic = [commonRoot, 'Input Data\VladislavSTD\Vladislav_STD\DataSet\depth']; %[commonRoot, 'Input Data/VladislavSTD/Vladislav_STD/depth'];   % 
    elseif dataSetNumber == 4
        depthPathBasic = [commonRoot, 'Input data\Mirela_dataset\Mirela_dataset'];
    elseif dataSetNumber == 5
        depthPathBasic = [commonRoot, 'Input data\Meshes\Aim@Shape_Selected'];
    end
    
    % before downsampling
    depthPath{1,1} = depthPathBasic;
    depthPath{2,1} = depthPathBasic;
    depthPath{3,1} = depthPathBasic;
    depthPath{4,1} = depthPathBasic;
%     depthPath{5,1} = [depthPathBasic, '_D1'];
%     depthPath{6,1} = [depthPathBasic, '_D1'];
%     depthPath{7,1} = [depthPathBasic, '_D2'];
%     depthPath{8,1} = [depthPathBasic, '_D2'];
    
    % after downsampling
    depthPath{1,2} = depthPathBasic;
    depthPath{2,2} = depthPathBasic;
    depthPath{3,2} = depthPathBasic;
%     depthPath{4,2} = [depthPathBasic, '_D1'];
%     depthPath{5,2} = [depthPathBasic, '_D1'];
%     depthPath{6,2} = [depthPathBasic, '_D2'];
%     depthPath{7,2} = [depthPathBasic, '_D2'];
%     depthPath{8,2} = [depthPathBasic, '_D3'];

    
end

