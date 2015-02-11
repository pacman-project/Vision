% this function is to ensure that learning ang recognition script use the
% same parametres

function [dxKernel, dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward, combs, largestLine, sigma, sigmaKernelSize, isErrosion, ...
                discRadius, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber) 

    if dataSetNumber == 1 || dataSetNumber == 3

        
        is_guided = true;
        r_guided = 4;
        eps = 0.05^2;
        is_mask_extended = false;
        maxExtThresh1 = 1;
        maxExtThresh2 = 60000;
        
        sigma = 0.8;
        sigmaKernelSize = 3;
        isErrosion = true;
        discRadius = 2;

    elseif dataSetNumber == 2

        % parameters for guided
        is_guided = false;
        r_guided = 10;
        eps = 0.05^2;
        % ----------------------
        is_mask_extended = true;
        maxExtThresh1 = 500;
        maxExtThresh2 = 900;
        sigma = 4; 
        sigmaKernelSize = round(2*sigma+1);
        isErrosion = true;
        discRadius = round(sigma - 1);

    end

    % parameters for prelimiminary processing functions
    dxKernel = load('dxKernel.mat');
    dxKernel = dxKernel.dxKernel;

    % helpful stuff for an inhobition function
    combs = load('settings/combs12.mat');
    combs = combs.combs; % combinations for the line discretization function
    largestLine = 12; % for decomposition of the line discretization function
    
    % forward and backward differences
    
    dyKernelTop = [-1, 0, 1, 0, 0]'/2;
    dyKernelBottom = [0, 0, -1, 0, 1]'/2;
    
    dxKernelBack  = [-1, 0, 1, 0, 0]/2;
    dxKernelForward = [0, 0, -1, 0, 1]/2;

end