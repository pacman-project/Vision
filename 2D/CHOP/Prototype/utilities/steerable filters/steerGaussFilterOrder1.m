function J = steerGaussFilterOrder1(I,theta,sigma,disFlag)
%    This function implments the steerable filter of the   
%    first derivative of Gaussian function of the paper
% 
%     W. T. Freeman and E. H. Adelson, "The Design
%     and Use of Steerable Filters", IEEE PAMI, 1991.
%
%    J = steerGaussFilterOrder1(I,theta,sigma,disFlag) 
%    calculated 
%    Input:
%    1. I: input image
%    2. theta: the orientation
%    3. sigma: standard deviation of the Gaussian template  
%    Output:
%    J. The response of derivative in theta direction
%
%    Author: Jincheng Pang, Tufts University, Dec. 2013.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part I: Assign algorithm parameters.


I = mean(double(I),3);

% Process input arguments (according to usage mode).
if ~exist('arg1','var') || ~isstruct(arg1)
   
    
    % Assign default filter orientation (if not provided).
    if ~exist('theta','var') || isempty(theta)
      theta = 0;
    end
    theta = -theta*(pi/180);

    % Assign default standard deviation (if not provided).
    if ~exist('sigma','var') || isempty(sigma)
       sigma = 1;
    end
    
    % Assign default visualization state (if not provided).
    if ~exist('disFlag','var') || isempty(disFlag)
       disFlag = false;
    end   
         
end % End of input pre-processing.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II: Evaluate separable filter kernels.

% Determine necessary filter support (for Gaussian).
Wx = floor((8/2)*sigma); 
if Wx < 1
  Wx = 1;
end
x = [-Wx:Wx];

[xx,yy] = meshgrid(x,x);
G1_0 = -(xx/sigma^2).*exp(-(xx.^2+yy.^2)/(2*sigma^2))/(sigma*sqrt(2*pi));
G1_90 = -(yy/sigma^2).*exp(-(xx.^2+yy.^2)/(2*sigma^2))/(sigma*sqrt(2*pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: Determine oriented filter response.

% Calculate image gradients (using separability).

Ix = imfilter(I,G1_0,'same','replicate');
Iy = imfilter(I,G1_90,'same','replicate');

% Evaluate oriented filter response.
J = cos(theta)*Ix+sin(theta)*Iy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part IV: Visualization

% Create figure window and display results.
% Note: Will only create output window is requested by user.
if disFlag
   figure(1); clf; set(gcf,'Name','Oriented Filtering');
   subplot(1,3,1); imagesc(I); axis image; colormap(gray);
      title('Input Image');
   subplot(1,3,2); imagesc(J); axis image; colormap(gray);
      title(['Filtered Image (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
   subplot(1,3,3); imagesc(cos(theta)*G1_0+sin(theta)*G1_90);
      axis image; colormap(gray);
      title(['Oriented Filter (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%