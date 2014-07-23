% this is to perform recognition of all levels of the hierarchy

% define folders configuration
root = 'D:\3D\LibHoP3D\';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'categorization']);
addpath([root,'Recognition']);

nClusters = 9;
combs = load('settings/combs10.mat');
combs = combs.combs; % combinations for the line discretization function
largestLine = 10; % for decomposition of the line discretization function

% ------------ parameters for line discretization--------------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
wCoverage = 0.25;  % parameters [0.01 - 0.5] make perfect sense 0.5 - more dense coverage
wOverlap = 0.1; % [0.01 - 0.3] 0.3 - supresses almost all overlaps

% [0.3, 0.3] - should give almost homogenious distribution

load('fx.mat');
strLen = length(fx);
output = zeros(1, strLen);

struct = load([root, 'settings/firstLayer']);
cluster1Centres = struct.cluster1Centres;
cluster1Lengths = struct.cluster1Lengths;
thresh = struct.thresh;

inds = 1:2:strLen-1;
fx = fx(inds);
strLen = length(fx);
[nearestClusters, errors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, cluster1Lengths, thresh); % replaces filter responces with nearest clusters
[outLine, ~] = lineDiscretizationOptLayer1(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine);

output(inds) = outLine;

kol = length(outLine(outLine>0))

% we have found out that for a regular case parameters 
% [wCoverage, woverlap] = [0.25, 0.05]

% for the areas where we need dense coverage [0.25, 0.001]
