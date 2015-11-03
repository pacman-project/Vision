%> Name: SetGaborParameters
%>
%> Description: This function sets the gabor parameters for feature
%> extraction.
%>
%> @param options Program options.
%>
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.08.2014
function [ options ] = SetGaborParameters( )
    options = [];
        %% ========== LOW - LEVEL FILTER PARAMETERS ==========
    options.numberOfGaborFilters = 8; % Number of Gabor filters at level 1.
    options.gaborFilterThr = 0.2; % Min response threshold for convolved features, 
                                  % taken as the percentage of max response 
                                  % in each image.
    options.absGaborFilterThr = 0; % Absolute response threshold for low-level 
                                   % responses.  
    options.gaborFilterSize = 10;       % Size of a gabor filter. Please note 
                                        % that the size also depends on the 
                                        % filter parameters, so consider them 
                                        % all when you change this!
    options.gabor.sigma = 1;            % Gabor filter parameters
    options.gabor.theta = 0;
    options.gabor.lambda = 1;
    options.gabor.psi = 0;
    options.gabor.gamma = 0.25;
    options.gabor.inhibitionRadius = 1;
                                        % The inhibition radius basically 
                                        % defines the half of the square's
                                        % size in which weaker responses other 
                                        % than the seed node will
                                        % be surpressed.
    options.gabor.stride = 1;           % Stride to use when extracting gabor
                                       % features.     
end

