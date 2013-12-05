% this is to perform recognition of all levels of the hierarchy

% define folders configuration
root = 'D:\3D\LibHoP3D\';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'categorization']);
addpath([root,'Recognition']);

% define the input data
depthPath = 'D:\3D\Input Data\Images for categorization\8view_1Scale\images';
depthPathDefault = [root, 'settings\list_depth.mat'];
fileListPrecomputed = true; 
is_subset = false; % whether we shall use all files for learning
subset_len = 100; % how much shall we use

% define and download the parameters
ss = load('settings/sigma.mat');
sigma = ss.sigma; % this is a sigma for gaussian derivative 
sigmaKernelSize = ss.sigmaKernelSize; % size of the kernel for gaussian
dxKernel = load('dxKernel.mat');
dxKernel = dxKernel.dxKernel;
nClusters = 9;
combs = load('settings/combs10.mat');
combs = combs.combs; % combinations for the line discretization function
largestLine = 10; % for decomposition of the line discretization function

% ------------ parameters for line discretization--------------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
wCoverage = 0.25;
wOverlap = 0.2; 

% --------Define the layers for recognition--------------------------------
is_first_layer = true;
is_second_layer = false; 
is_third_layer = false;
is_4th_layer = false; 
is_5th_layer = false; 
is_6th_layer = false;

areLayersRequired = [is_first_layer, is_second_layer, is_third_layer, is_4th_layer, is_5th_layer, is_6th_layer];
% --------output files-----------------------------------------------------

outRoot = 'D:\3D\Input Data\Images for categorization\1view_1Scale\ElementsOpt';

%--------------------------------------------------------------------------
% forming a filelist here

[list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
struct = load([root, 'settings/firstLayer']);
cluster1Centres = struct.cluster1Centres;
cluster1Lengths = struct.cluster1Lengths;
thresh = struct.thresh;

%--------------------------------------------------------------------------

computeCoverage(list_depth, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, cluster1Centres, cluster1Lengths, thresh, ...
                areLayersRequired, outRoot, combs, largestLine, wCoverage, wOverlap);





