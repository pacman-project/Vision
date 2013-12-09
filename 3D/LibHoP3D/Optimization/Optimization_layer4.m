% this is to compute the layer 3 with optimization task

% optimization procedure minimizes the following equation:
% given set of parts S (triples)
% select subset S' such that
% minimize [ sum(S - S'(nearest)) + alpha*betweenDist  + betta*length(S') + gamma*skewedness(S') ]

% where betweenDist - distances between S'

function [outTriple] = Optimization_layer4(X, frequencies, triples3Out, triples4, nClusters, n3Clusters, numDisp, subset_len)

    disp('Optimization of the 4th layer started...');

    adopted = [];
    alpha = 0.2 * sum(frequencies) / subset_len;
    betta = sum(frequencies) / subset_len;  % to penalize number of clusters  2000
    gamma = sum(frequencies) / subset_len;  % to penalize curvedness (and skewedness) of the clusters  1000
    
    distThresh = 2;

    for i = 1 : nClusters
        for j = 1 : nClusters
            cur = [i, j, i, j, i, j, i, j, i, j, i, j, i, j, i, j, i, j]; % flat elements of all different orientations
            adopted = [adopted; cur];
        end
    end

    len = length(adopted);
    
    % here we implement random walk from each adopted elements
    added = [];

    for j = 0 : 40

        for k = 1:2
            if (j == 0 && k == 1)
                i = 41;
            elseif (j == 0 && k == 2)
                continue;
            else
                if k == 1
                    i = 41 - j;
                end
                if k == 2
                    i = 41 + j;
                end
            end

            str = ['      ', num2str(adopted(i,:))];
            disp(str);


            while true
                traversal_done = false;
                curEl = adopted(i,:);
                score = computeScore(X, frequencies, [adopted; added], alpha, betta, distThresh);    
                [elNext, scoreOut, is_found] = steepestSlope4(X, triples4, frequencies, curEl, score, triples3Out, n3Clusters, adopted, added, alpha, betta, gamma, numDisp, distThresh);
                if is_found == 0
                    break;
                end

                % traverse to this direction until minimum is reached
                while is_found
                    traversal_done = true;
                    elReserve = elNext;
                    [elNext, scoreOut, is_found] = steepestSlope4(X, triples4, frequencies, elNext, scoreOut, triples3Out, n3Clusters, adopted, added, alpha, betta, gamma, numDisp, distThresh);
                end

                % traversal ended
                added = [added; elReserve];

                str = ['added ', num2str(elReserve)];
                disp(str);
            end
        end
    end
    
    outTriple = [adopted; added];
end
        




