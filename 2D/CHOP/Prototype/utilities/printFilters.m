%> Name: printFilters
%>
%> Description: Given square filters in a 
%>
%> @param W NxD array which contains N filters.
%> @param dimensions 1xM array, which contains the filter dimensions.
%> M is either 2 or 3.
%> @param datasetName The name of the dataset which this filters are
%> trained for. 
%>  
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.08.2015
function [] = printFilters( W, dimensions, datasetName )
    filterFolderName = [pwd '/filters/auto'];
    
    for filterItr = 1:size(W,1)
            gaborFilt = reshape(W(filterItr,:), dimensions);                                             
            printedFilt = uint8(round(255 * (gaborFilt - min(min(min(gaborFilt)))) / (max(max(max(gaborFilt))) - min(min(min(gaborFilt))))));
            imwrite(printedFilt, [filterFolderName '/filt' num2str(filterItr) '.png']);
            save([filterFolderName  '/filt' num2str(filterItr) '.mat'], 'gaborFilt');
    end
    
    C = [];
    mu = [];
    invMat = [];
    whMat = [];
    if ~exist([pwd  '/filters/vis/' datasetName], 'dir')
         mkdir([pwd  '/filters/vis/' datasetName]);
    end
    save([pwd  '/filters/vis/' datasetName '/C.mat'], 'C', 'mu', 'invMat', 'whMat');
end

