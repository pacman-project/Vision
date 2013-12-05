% this is the main function to learn statistics of the 3th layer

clear all;

sigma = 0.6;
kernelSize = 3;
thresh = 95;  % for perspective camera
halfThresh = thresh/2;
depthStep = thresh/3;
nClusters = 9;
n3Clusters = nClusters^2;
smallestLine = 20; % the smallest line used for learning
subset_len = 6000;


fieldSize = 19;
halfFieldSize = floor(fieldSize/2);
depthSamples = 71;  %shows depths differences [-25*depthStep .. 25*depthStep]
halfDepthSamples = floor(depthSamples / 2); % should be 25

windowMarks = zeros(fieldSize, fieldSize);
marks1 = zeros(fieldSize, fieldSize);

windowI = zeros(fieldSize,fieldSize);
windowIx = zeros(fieldSize,fieldSize);
windowIy = zeros(fieldSize,fieldSize);

dxKernel4 = [1,-8,0,8,-1]; % consistancy order 4 
dxKernel4 = dxKernel4 / 12;
shiftX = 0;
shiftY = 0;

kernelW = load('settings/kernelW.mat');
kernelW = kernelW.kernelW;    % this is the kernel for lineDiscretization and lineDiscretizationY

% cluster sizes in percent (if we have nClusters = 9)
if nClusters == 9
    clusterSizes = load('settings/cluster1Sizes.mat'); 
    clusterSizes = clusterSizes.clusterSizes;
    cluster1SizesPercent = clusterSizes* 0.01;
end

[clusterCenters, cluster1Length] = defineCluster1Centers(nClusters, cluster1SizesPercent, thresh);

strOutput = 'D:/3D/Elements 2Layer/';
strOutE = 'Elements/';
strOutTI = 'TrimedImages/';
strRootE = [strOutput, strOutE];
strRootTI = [strOutput, strOutTI];
fileListPrecomputed = true;
is_subset = true; % if we take only a subset of images

statistics = zeros(n3Clusters, n3Clusters, fieldSize, fieldSize, depthSamples);

if fileListPrecomputed == true   % this is done for speed up reasons
    list_el = load('listEl.mat');
    list_el = list_el.list_el;
    list_TR = load('listTR.mat');
    list_TR = list_TR.list_TR;
else
    list_el = load_filelist(strRootE);
    list_TR = load_filelist(strRootTI);
end 

lenF = length(list_el);
if (is_subset == true)
    nums = randperm(lenF);
    nums = nums(1:subset_len); 
    lenF = subset_len;
end

for i = 1:lenF % To use later for models
        if (is_subset == true)
            marks = imread(list_el{nums(i)});
            IG = imread(list_TR{nums(i)});
        else
            marks = imread(list_el{i});
            IG = imread(list_TR{i}); % list should contain the full path
        end
        IG = double(IG);        
        [r,c] = size(IG);
        
        % create a mask
        [r,c] = size(IG);
        mask = zeros(r,c);
        mask(IG > 0) = 1;
        
        % check if the mask is empty or not
        maxM = max(max(mask));
        if maxM == 0 
            continue;
        end   
        
%         [I, mask] = trimImage(I, mask); % trim the image
        
        % compute filter responces
        hx = dxKernel4;
        hy = hx';       
        
%         h = fspecial('gaussian', kernelSize, sigma);
%         IG = imfilter(I, h); % gaussian                 % done

        % now computing Gaussian derivatives
        Ix = imfilter(IG, hx);
        Ix = Ix.*mask; % to avoid high on the boundary
        Iy = imfilter(IG, hy); % derivative in y direction
        Iy = Iy.*mask;
        
%         % perform line discretization (line after line)
%         [r, c] = size(Ix);
%         marks = zeros(r,c); % here we collect marks
%         for k = 1:r % for every line we perform discretization
%             ind = find(mask(k,:));
%             if length(ind) < smallestLine
%                 marks(k,:) = zeros(1, c); % line of zeros in this case
%             else
%                 fx = Ix(k,ind(1):ind(end));
%                 strlen = length(fx); % strlen must be greater than smallestLine
%                 output = lineDiscretization(fx, strlen, nClusters, clusterCenters, cluster1Length, thresh);
%                 outline = zeros(1,c);
%                 marks(k,ind(1):ind(end)) = output;
%             end
%         end
%        % imtool(marks, [1,nClusters]);  %   GIVES NICE PICTURES 
%         marks = marks .* mask;
        
        marks = addZeroBoundaries(marks, halfFieldSize);      
        I = addZeroBoundaries(IG, halfFieldSize);
        Ix = addZeroBoundaries(Ix, halfFieldSize);
        Iy = addZeroBoundaries(Iy, halfFieldSize);
        mask = addZeroBoundaries(mask, halfFieldSize);
        
        [rows,cols] = find(marks>0);
        ll = length(rows);
        
        % we now put every element to the central position of the windows
        for  j = 1:ll
            
            if marks(rows(j), cols(j)) == 41
                a = 2;
            end
            
            windowI =  I(rows(j)-halfFieldSize : rows(j)+halfFieldSize, cols(j)-halfFieldSize : cols(j)+halfFieldSize);
            windowIx = Ix(rows(j)-halfFieldSize : rows(j)+halfFieldSize, cols(j)-halfFieldSize : cols(j)+halfFieldSize);
            windowIy = Iy(rows(j)-halfFieldSize : rows(j)+halfFieldSize, cols(j)-halfFieldSize : cols(j)+halfFieldSize);
            windowMask = mask(rows(j)-halfFieldSize : rows(j)+halfFieldSize, cols(j)-halfFieldSize : cols(j)+halfFieldSize);
            
            % discretize every row
            for jj = 1:fieldSize
                output = lineDiscretization(windowIx(jj,:), fieldSize, nClusters, kernelW, clusterCenters, cluster1Length, thresh);
                windowMarks(jj,:) = output;
            end
            
            windowMarks = windowMarks.*windowMask; 
            marks1 = zeros(fieldSize, fieldSize);

            
            % next we create second layer elements in every window
            % we predict centers in the positions: halfFieldSize + 1,
            % halfFieldSize + 1 - 5, halfFieldSize + 1 + 5
            
            posY = [halfFieldSize - 4, halfFieldSize + 1, halfFieldSize + 6];
            width = 2;
  
            for kk = 1:3
                broadLine = windowMarks(:, posY(kk) - 2:posY(kk) + 2);
                narrowLine = zeros(fieldSize, 1);
                
                for iii = 1:fieldSize     % put everything from broad line to narrow one
                    rr = broadLine(iii,:);
                    inds = find(rr>0);
                    if ~isempty(inds)
                        narrowLine(iii) = broadLine(iii, inds(1));
                    end
                end
                
                if length(find(narrowLine>0)) >= 5
                    outline = lineDiscretizationY(narrowLine, fieldSize, kernelW, nClusters);
                
                    inds = find(outline > 0);
                    inds = inds(inds > 2);
                    inds = inds(inds < 18);
                    for iii = 1:length(inds)
                        indC = [inds(iii) , posY(kk)];
                        central = outline(inds(iii));                   
                        
                        window = windowMarks(indC(1)-2:indC(1)+2, indC(2)-2:indC(2)+2);               
                        % adjust the element in x direction
                        [rr,cc] = find(window == central | window == central-1 | window == central+1);
                        shiftX = - 3 + round(sum(cc)/length(cc));

                        dY = windowIy(indC(1), indC(2) + shiftX);
                        clusterY = define1Cluster(dY, nClusters, cluster1Length, thresh);
                        clusterXY = compute2elementIndex(central, clusterY, nClusters);
                        marks1(indC(1), indC(2) + shiftX) = clusterXY;
                    end
                else
                    continue;
                end
            end
%                 
%           % collect the statistics!
            % adjust central element to the center of the field of view
            center = halfFieldSize+1;
            
            wind = marks1(center-1:center+1, center-1:center+1);
            [rw,cw] = find(wind>0);
            if length(rw) > 1
                % something is fully wrong
                continue;
            elseif length(rw) == 1
                shiftX = cw(1) - 2;
                shiftY = rw(1) - 2;
            elseif length(rw) == 0
                % then take a larger window
                wind = marks1(center-2:center+2, center-2:center+2);
                [rw,cw] = find(wind>0);                    
                if length(rw) ~= 1 % still
                    continue;
                else
                    shiftX = cw(1) - 3;
                    shiftY = rw(1) - 3;
                end
            end

            % now we pretend we know shiftX and shiftY
            el = marks1(halfFieldSize + 1 + shiftY, halfFieldSize + 1 + shiftX);
            depthCentral = windowI(halfFieldSize + 1 + shiftY, halfFieldSize + 1 + shiftX);
            marks1(halfFieldSize + 1 + shiftY, halfFieldSize + 1 + shiftX) = 0;
            

            % find all the elements in the window (except the central)
            [rw,cw] = find(marks1>0);
            if isempty(rw)
                continue;
            end
            for nn = 1:length(rw)         
                curEl = marks1(rw(nn), cw(nn));
                curDepth = windowI(rw(nn), cw(nn));
                shiftDepth = curDepth - depthCentral;
                
                % now quantize this value
                shD = int32(round(shiftDepth / depthStep));
                if shD > halfDepthSamples
                    shD = halfDepthSamples;
                elseif shD < - halfDepthSamples
                    shD = - halfDepthSamples;
                end
                shDp = halfDepthSamples + 1 + shD; % positive
                shDn = halfDepthSamples + 1 - shD; % negative
                
                statistics(el, curEl, cw(nn) + shiftX, rw(nn) + shiftY, shDp) = statistics(el, curEl, cw(nn) + shiftX, rw(nn) + shiftY, shDp) + 1;
%               statistics(curEl, el, cw(nn) + shiftX, rw(nn) + shiftY, shDn) = statistics(curEl, el, cw(nn) + shiftX, rw(nn) + shiftY, shDn) + 1;
 
            end         
            a = 1;        
        end
          
%     if mod(i,10) == 0
%         i
%     end
    i
end

save('statistics/statistics5D_3layer_19_6000.mat','statistics');












