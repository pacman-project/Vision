% this is to compute the layer 3 with optimization task

% optimization procedure minimizes the following equation:
% given set of parts S (triples)
% select subset S' such that
% minimize [ sum(S - S'(nearest)) + alpha*betweenDist  + betta*length(S') + gamma*skewedness(S') ]

% where betweenDist - distances between S'

function [triples4Out, n4Clusters] = Optimization_layer4(X, frequencies, triples3Out, triples4, nClusters, n3Clusters, maxDist, ...
                                                            subset_len, alphaParam, betaParam, gammaParam)

    disp('Optimization of the 4th layer started...');

    adopted = [];
    alpha = alphaParam * sum(frequencies) / subset_len;
    betta = betaParam  * sum(frequencies) / subset_len;  % to penalize number of clusters  
    gamma = gammaParam * sum(frequencies) / subset_len;  % to penalize curvedness (and skewedness) of the clusters 
    
    distThresh = 2;

%     for i = 1 : nClusters
%         for j = 1 : nClusters
%             cur = [i, j, i, j, i, j, i, j, i, j, i, j, i, j, i, j, i, j]; % flat elements of all different orientations
%             adopted = [adopted; cur];
%         end
%     end

    for i = 1:n3Clusters
        cur = triples3Out(i,:);
        cur = [cur, cur, cur];
        adopted = [adopted; cur];
    end

    len = length(adopted);
    
    % in order to speed up all things we have to pre-compute all distances
    % between layer 3 parts
    rx = size(X, 1);
    
    distXtoInitAdopted = Isodata_distances(X, adopted, rx, len, false, false); % distances form X to initially adoped elements
    
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

            counter = 0;
            while true
                traversal_done = false;
                curEl = adopted(i,:);
                score = computeScore(X, frequencies, adopted, added, alpha, betta, distThresh, distXtoInitAdopted);    
                [elNext, scoreOut, is_found] = steepestSlope4(X, triples4, frequencies, curEl, score, triples3Out, n3Clusters, ...
                                            adopted, added, alpha, betta, gamma, maxDist, distThresh, distXtoInitAdopted);
                if is_found == 0
                    break;
                end

                % traverse to this direction until minimum is reached
                while is_found
                    traversal_done = true;
                    elReserve = elNext;
                    [elNext, scoreOut, is_found] = steepestSlope4(X, triples4, frequencies, elNext, scoreOut, triples3Out, n3Clusters, ...
                                           adopted, added, alpha, betta, gamma, maxDist, distThresh, distXtoInitAdopted);
                end

                % traversal ended
                added = [added; elReserve];
                counter = counter + 1;
                str = ['added ', num2str(elReserve)];
                disp(str);
                str = ['size of vocabulary is now ',num2str(size([adopted; added], 1))];
                disp(str);
                
                if counter >= 5  % merge added and adopted
                    lenAdded = size(added, 1);
                    adopted = [adopted; added];
                    distAdder = Isodata_distances(X, added, rx, lenAdded, false, false);
                    distXtoInitAdopted = [distXtoInitAdopted, distAdder];
                    added = [];
                    counter = 0;
                end

            end
        end
    end
    
    outTriple = [adopted; added];
    n4Clusters = size(outTriple, 1);
    triples4Out = zeros(n4Clusters, 3);
    
    table3 = uint16(zeros(9,9,9,9,9,9));
    
    for i = 1:n3Clusters
        cur = triples3Out(i,:);
        table3(cur(1),cur(2),cur(3),cur(4),cur(5),cur(6)) = i;
    end
    
    % convert them to format top, middle, bottom
    for i = 1:n4Clusters
        for j = 0:2
            cur = outTriple(i, j*6+1 : (j+1)*6);
            % find the corresponding 3rd layer element
            triples4Out(i,j+1) = table3(cur(1),cur(2),cur(3),cur(4),cur(5),cur(6));
        end
    end
   
end
        





