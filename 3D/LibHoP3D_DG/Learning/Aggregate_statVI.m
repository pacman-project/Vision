% this is to aggregate statistics for the third layer and prepare it for
% optimization function
% triples should look like that: left, middle, right

function [X, frequencies, triples] = Aggregate_statVI(lenFiles, nPrevClusters, sieve_thresh, pairsLeft, pairsRight) 

    is_visualization = false;
   
    triples = zeros(nPrevClusters, nPrevClusters, nPrevClusters);
    disp('build triples...');
    curTS = 200;
    
    for k = 1:lenFiles %20:20:lenFiles
        
        strIn = ['Temp/outList_',num2str(k), '.mat'];
        if ~exist(strIn, 'file')
            continue;
        end
        aa = load(strIn);
        outList = aa.outList;
        
        ids1 = outList(1, :) > 0;
        ids3 = outList(3, :) > 0;
        ids = ids1 & ids3;
        
        outList = outList(:, ids);
        
        for i = 1:size(outList, 2)
            triples(outList(1,i), outList(2,i), outList(3,i)) = triples(outList(1,i), outList(2,i), outList(3,i)) + 1;
        end
    end

    disp('aggregate statistics');
    ind = 0;
    X = zeros(curTS,3);
    frequencies = zeros(curTS,1);

    for i = 1:nPrevClusters
        for j = 1:nPrevClusters
            for k = 1:nPrevClusters
%                 if ~is_sparse
                    if triples(i,j,k) >= sieve_thresh
                        ind = ind + 1;

                        X(ind,1) = i;  % left   or top
                        X(ind,2) = j;  % centre  or centre
                        X(ind,3) = k;  % right   or bottom
                        frequencies(ind) = triples(i,j,k); 
                    end
%                 else  % if is_sparse
%                     if triples{i}(j,k) >= sieve_thresh
%                         ind = ind + 1;
% 
%                         X(ind,1) = i;  % left   or top
%                         X(ind,2) = j;  % centre  or centre
%                         X(ind,3) = k;  % right   or bottom
%                         frequencies(ind) = triples{i}(j,k); 
%                     end
%                 end
            end
        end
    end
    X = X(1:ind, :);
    frequencies = frequencies(1:ind);
    
    [frequencies, iis] = sort(frequencies, 'descend');
    X= X(iis, :);

    
    % visualization of the results
    if is_visualization
        for j = 1:ind

            figure;
            disp(frequencies(j));

            quiver3(0,0,0, 0,0,0.01);
            xlim([-0.02 0.02])
            ylim([-0.02 0.02])
            zlim([-0.02 0.02])
            hold on

            for i = 1:2
                if i == 1
                    partsV = pairsLeft(X(j, 1), 3:end);
                    partsV(4:6) = partsV(4:6)/100;
                elseif i == 2
                    id = X(j, 3) - size(pairsLeft, 1);
                    partsV = pairsRight(id ,3:end);
                    partsV(4:6) = partsV(4:6)/100;
                end
                quiver3(partsV(1), partsV(2), partsV(3), partsV(4), partsV(5), partsV(6), 'color', 'blue');
                hold on
            end

            a = 2;
        end

        allPlots = findall(0, 'Type', 'figure');
        delete(allPlots);
    end
    
    strOut = 'Temp/Aggregated.mat';
    save(strOut, 'X', 'frequencies', 'triples');
end

