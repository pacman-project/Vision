function [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, optimizationFlag)

    optimizationFilterSizes = [21 21 21 21 21 31 41 61 71];

    dummySigma = 0.15;
    if optimizationFlag
         angleStep = 3; % 22.5 for no rotations!
    else
         angleStep = 180/numel(options.filters);
    end
    
    minPixelValue = 1/255;
    
    % Obtain visual filter.
    firstVisFilter = options.filters{1};
    firstVisFilter = (firstVisFilter - min(min(firstVisFilter))) / (max(max(firstVisFilter)) - min(min(firstVisFilter)));
    
    % Generate rotated filters.
    optimizedFilters = cell(numel(optimizationFilterSizes),1);
    numberOfOptFilters = round(180/angleStep);
    
    % Generate barrel filter (initial)
    padSize = 5;
    for optFiltSizeItr = 1:numel(optimizationFilterSizes)
         filterSize = optimizationFilterSizes(optFiltSizeItr) + 2*padSize;
         vals = normpdf((0:(filterSize-1))/(filterSize-1), 0.5, 0.125);
         vals = (vals-min(vals))/(max(vals) - min(vals));
         firstFilter = repmat(vals, filterSize, 1);
         curOptimizedFilters = cell(numberOfOptFilters,1);
         for filterItr = 0:((180-angleStep)/angleStep)
              curAngle = -angleStep * filterItr;
              curFilter = imrotate(firstFilter, curAngle, 'bilinear', 'crop');
              curFilter = curFilter((padSize+1):(end-padSize), (padSize+1):(end-padSize));
              curFilter(curFilter<minPixelValue) = minPixelValue;
              curFilter = round(curFilter * 255)/255;
              curOptimizedFilters{filterItr+1} = curFilter;
         end
         
         %% Now, we convert the filters to uint8.
         for filterItr = 1:numel(curOptimizedFilters)
              curOptimizedFilters{filterItr} = uint8(round(curOptimizedFilters{filterItr} * 255));
         end
         optimizedFilters{optFiltSizeItr} = double(cat(3, curOptimizedFilters{:}));
    end
    
    %% Create filters for visualization.
    visFilters = cell(numberOfOptFilters,1);
    for filterItr = 0:((180-angleStep)/angleStep)
         curAngle = -angleStep * filterItr;
         curVisFilter = imrotate(firstVisFilter, curAngle, 'bilinear', 'crop');
         curVisFilter(curVisFilter<minPixelValue) = minPixelValue;
         curVisFilter = round(curVisFilter * 255)/255;
         visFilters{filterItr+1} = curVisFilter;
    end
        
    %% Now, we convert the filters to uint8.
    for filterItr = 1:numel(visFilters)
         visFilters{filterItr} = uint8(round(visFilters{filterItr} * 255));
    end
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
    rfSizes = zeros(10,1);
    for rfSizeItr = 1:min((size(rfSizes,1)-1), numel(options.receptiveFieldSizes))
         rfSize = getRFSize(options, rfSizeItr+1);
         rfSizes(rfSizeItr+1) = rfSize(1);
    end
end