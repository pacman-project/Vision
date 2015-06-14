% inference based on co-occurrence statistics

function [] = InferenceNext_simple(statistics, outputCoords, outputScales, outputFrames, curTS, inFolder, outFolder, list_El, triplesCurOut, NPrevClusters, nNClusters)

    lenF = length(list_El);
    firstColumn = outputCoords(:,1);
    triplesCurOut = triplesCurOut(1:nNClusters, :);
    
    iDDs = [2,1,4];
    statistics = statistics(:, iDDs);
    lenDPW = length(inFolder);
    
    triples = zeros(NPrevClusters,NPrevClusters,NPrevClusters);
    
    for i = 1:nNClusters
        triples(triplesCurOut(i,1), triplesCurOut(i,2), triplesCurOut(i,3)) = i;
    end
    
    for i = 1:lenF
         curStr = list_El{i};
         marksPrev = imread(list_El{i});
         [r,c] = size(marksPrev);
         
         marksOut = zeros(r,c);
         
         % Select all parts which refer to this image
         ind = find(firstColumn == i);
         lenE = length(ind);
         
         if lenE == 0
             continue;
         end
         
         curStat = statistics(ind, :);
         curCoords = outputCoords(ind, :);
         lenS = size(curStat, 1);
         
         
        for j = 1:lenS

            curTriple = curStat(j, :);
            partID = triples(curTriple(1), curTriple(2), curTriple(3));

            if partID ~= 0
             x = curCoords(j,2);
             y = curCoords(j,3);
             marksOut(y,x) = partID;
            end
        end
        
        % write the result to a file 
        fileName = curStr(lenDPW+1:end);
        outFile = [outFolder, fileName];

        ll = strfind(outFile, '/');
        lll = strfind(outFile, '\');
        ll = [ll, lll];
        ll = max(ll); % last position
        folderName = outFile(1:ll);


        b = exist(folderName,'dir');
        if b == 0
            mkdir(folderName);
        end
        b = exist(outFile, 'file');

        marksOut = uint8(marksOut);  
    end

end

