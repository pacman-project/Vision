function [] = addPaths(root)

    addpath([root,'utility']);
    addpath([root,'utility/GuidedFilter']);
    addpath([root,'settings']);
    addpath([root,'Optimization']);
    addpath([root,'Learning']);
    addpath([root,'visualization']);
    addpath([root,'statistics']);
    addpath([root,'Temp']);
    addpath([root,'ISODATA']);
    addpath([root,'Inference']);
    addpath([root,'categorization']);
    addpath([root,'DifferentialGeometry']);
    addpath([root,'toolbox_graph\toolbox_graph']);
    addpath([root,'toolbox_graph\toolbox_graph\toolbox']);
    addpath([root,'quaternions-1.3\quaternions'])
    addpath([root,'quaternions-1.3\quaternions\test'])
    addpath([root,'toolbox_signal']);
    addpath([root,'toolbox_general']);
    addpath([root,'toolbox_wavelet_meshes']);
    addpath([root,'toolbox_fast_marching']);
    addpath([root,'CircStat2012a']);
end

