% this function is to ensure that learning ang recognition script use the
% same parametres

function [dxKernel, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber) 

    if dataSetNumber == 1 || dataSetNumber == 3

        sigma = 1.2;
        sigmaKernelSize = 5;
        isErrosion = true;
        discRadius = 2;
        
        is_guided = false;
        r_guided = 1;
        eps = 0.1;
        is_mask_extended = false;
        maxExtThresh1 = 1;
        maxExtThresh2 = 60000;

    elseif dataSetNumber == 2

        % parameters for guided
        is_guided = false;
        r_guided = 10;
        eps = 0.05^2;
        % ----------------------
        is_mask_extended = true;
        maxExtThresh1 = 500;
        maxExtThresh2 = 900;

        % parameters for gaussian
        sigma = 5;
        sigmaKernelSize = round(2*sigma+1);
        % -----------------------

        isErrosion = true;
        discRadius = 6;
    end

    % parameters for prelimiminary processing functions
    dxKernel = load('dxKernel.mat');
    dxKernel = dxKernel.dxKernel;

    % helpful stuff for an inhobition function
    combs = load('settings/combs12.mat');
    combs = combs.combs; % combinations for the line discretization function
    largestLine = 12; % for decomposition of the line discretization function

end