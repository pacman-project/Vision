% this function is to ensure that learning ang recognition script use the
% same parametres

function [filtOptions, zeroThresh]  = loadFilteringParameters(dataSetNumber) 

    if dataSetNumber == 1 || dataSetNumber == 3

        filtOptions.is_guided = true;
        filtOptions.r_guided = 8;
        filtOptions.eps = 0.3^2;
        filtOptions.is_mask_extended = false;
        filtOptions.minDepth = 1;
        filtOptions.maxDepth = 60000;
        
%         filtOptions.sigma = 0.8;
%         filtOptions.sigmaKernelSize = 3;
        filtOptions.isErrosion = true;
        filtOptions.discSize = 2;
        
        zeroThresh = 10^(-3);

    elseif dataSetNumber == 2

        % parameters for guided
        filtOptions.is_guided = false;
%         filtOptions.r_guided = 10;
%         filtOptions.eps = 0.05^2;
        % ----------------------
        filtOptions.is_mask_extended = true;
        filtOptions.minDepth = 500;
        filtOptions.maxDepth = 900;
        filtOptions.sigma = 4; 
        filtOptions.sigmaKernelSize = round(2*sigma+1);
        filtOptions.isErrosion = true;
        filtOptions.discRadius = round(sigma - 1);

        
    elseif dataSetNumber == 4

        % parameters for guided
        filtOptions.is_guided = false;
        filtOptions.r_guided = 10;
        filtOptions.eps = 0.05^2;
        % ----------------------
        filtOptions.is_mask_extended = true;
        filtOptions.minDepth = 500;
        filtOptions.maxDepth = 1100;
        filtOptions.sigma = 3.5; 
        filtOptions.sigmaKernelSize = round(2*sigma+1);
        filtOptions.isErrosion = true;
        filtOptions.discRadius = round(sigma - 1);

    end

    % parameters for prelimiminary processing functions
    dxKernel = load('dxKernel.mat');
    filtOptions.dxKernel = dxKernel.dxKernel;
    
    % forward and backward differences 
    filtOptions.dyKernelTop = [-1, 0, 1, 0, 0]'/2;
    filtOptions.dyKernelBottom = [0, 0, -1, 0, 1]'/2;
    filtOptions.dxKernelBack  = [-1, 0, 1, 0, 0]/2;
    filtOptions.dxKernelForward = [0, 0, -1, 0, 1]/2;
     
    filtOptions.isX = true;
    filtOptions.isY = true;
    filtOptions.isX_FB = false;
    filtOptions.isY_FB = false;
    filtOptions.dilateMaskDeriatives = true;
    filtOptions.deriativesThresh = 6.0;
    filtOptions.isTrim = false;
    filtOptions.isFiltering = true;
    
end