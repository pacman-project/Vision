% this function is to visualize any 4th layer elements. It basically calls
% the function surfaceVisualizer(fieldSize, positions, elements, nClusters, thresh, depthStep)

% positions: [x1,y1,z1; x2,y2,z2; ...] may be real numbers far from margin
% elements (35,44,51,...)
% nCluster - integer
% thresh - real number (refers to thresh form the 1st layer)
% fieldSize is a vector [sizeX, SizeY, sizeZ]
% depthStep - real number (presumably 95/3)


% elements must contain 9 numbers ex. [41,41,41, 41,41,41, 41,41,41]


function [ out ] = VisualizerLayer4( elements, nClusters, fieldSize, depthStep, thresh)

    stats = load('statistics/pairs3layer_8.mat');
    stats = stats.cluster3stats;
    
    if (nargin < 5)
        thresh = 95;
    end
    if (nargin < 4)
        depthStep = thresh / 3;
    end
    if (nargin < 3)
        fieldSize = [19, 19, 71];
    end
    if (nargin < 2)
        nClusters = 9;
    end
    if (nargin < 1)
        elements = [41,41,41, 41,41,41, 41,41,41];
    end
    
    % we have to define positions of the elements here
    line1 = stats(:,1);
    line2 = stats(:,2);
    line16 = stats(:,16);
    
    positions = zeros(9,3);
    positions(1,:) = [10, 10, 36];
    
    for i = 1:8
        disp = i + 1;
        lineIndex = find(line1 == elements(1) & line2 == elements(disp) & line16 == disp);
        
        if isempty (lineIndex)
            elements(disp) = [];
            continue;
        end
        line = stats(lineIndex,:);
        mu = line(3:5);
        mu = round(mu);
        positions(disp,:) = mu;
    end
   
    if ~is_failed
        surfaceVisualizer(fieldSize, positions, elements, nClusters, thresh, depthStep);
    end
    
    out = ~is_failed;
end

