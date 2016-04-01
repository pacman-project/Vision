function [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, optimizationFlag)

    dummySigma = 0.15;
    if optimizationFlag
         angleStep = 3; % 22.5 for no rotations!
    else
         angleStep = 180/numel(options.filters);
    end
    
    if strcmp(options.filterType, 'gabor')
        stride = options.gabor.stride;
    else
        stride = options.auto.stride;
    end
    
    minPixelValue = 1/255;
    
    % As a case study, we replace gabors with 1D gaussian filters
    % stretched.
    visFilterSize = size(options.filters{1},1);
    
    % Generate barrel filter (initial)
    filterSize = visFilterSize + 16;
    padSize = 5;
    vals = normpdf(1:filterSize, (filterSize+1)/2, 5);
    vals = vals/max(vals);
    firstFilter = repmat(vals, filterSize, 1);
    
    % Obtain visual filter.
    firstVisFilter = options.filters{1};
    firstVisFilter = (firstVisFilter - min(min(firstVisFilter))) / (max(max(firstVisFilter)) - min(min(firstVisFilter)));
    
    % Generate rotated filters.
    numberOfOptFilters = round(180/angleStep);
    optimizedFilters = cell(numberOfOptFilters,1);
    visFilters = cell(numberOfOptFilters,1);
    for filterItr = 0:((180-angleStep)/angleStep)
         curAngle = -angleStep * filterItr;
         curFilter = imrotate(firstFilter, curAngle, 'bilinear', 'crop');
         curFilter = curFilter((padSize+1):(end-padSize), (padSize+1):(end-padSize));
         curVisFilter = imrotate(firstVisFilter, curAngle, 'bilinear', 'crop');
         curFilter(curFilter<minPixelValue) = minPixelValue;
         curVisFilter(curVisFilter<minPixelValue) = minPixelValue;
         curFilter = round(curFilter * 255)/255;
         curVisFilter = round(curVisFilter * 255)/255;
         optimizedFilters{filterItr+1} = curFilter;
         visFilters{filterItr+1} = curVisFilter;
    end
        
    %% Now, we convert the filters to uint8.
    for filterItr = 1:numel(optimizedFilters)
         visFilters{filterItr} = uint8(round(visFilters{filterItr} * 255));
         optimizedFilters{filterItr} = uint8(round(optimizedFilters{filterItr} * 255));
    end
    optimizedFilters = double(cat(3, optimizedFilters{:}));
    visFilters = double(cat(3, visFilters{:}));
    
    %% Here, we create a likelihood lookup table for predictions.
    likelihoodLookupTable = zeros(256);
    startVals = (-1/510:1/255:(1+1/510))';
    for itr = 0:255
         probs = normcdf(startVals, ones(257, 1) * (itr/255) , repmat(dummySigma,257,1));
 %        normProb = probs(end) - probs(1);
         normProb = 1;
         likelihoodLookupTable(itr+1,:) = (probs(2:end) - probs(1:(end-1))) / normProb;
    end
    likelihoodLookupTable = log(likelihoodLookupTable);
    rfSizes = zeros(levelItr,1);
    for rfSizeItr = 2:size(rfSizes,1)
         rfSize = getRFSize(options, rfSizeItr);
         rfSizes(rfSizeItr) = rfSize(1);
    end
end