% this is to compute the layer 3 with optimization task

% optimization procedure minimizes the following equation:
% given set of parts S (triples)
% select subset S' such that
% minimize [ sum(S - S'(nearest)) + alpha*betweenDist  + betta*length(S') + gamma*skewedness(S') ]

% where betweenDist - distances between S'

function [outTriple] = Optimization_layer3(X, frequencies, triples, nClusters)

    disp('Optimization started...');

    adopted = [];
    alpha = 0;
    betta = sum(frequencies) / 3000;  % to penalize number of clusters  2000
    gamma = sum(frequencies) / 3000;  % to penalize curvedness (and skewedness) of the clusters  1000
    distThresh = 1;

    for i = 1 : nClusters
        for j = 1 : nClusters
            adopted = [adopted; i, j, i, j, i, j];
        end
    end

    len = length(adopted);

    displacements = [
        -1, 0, 0, 0, 0, 0;
         1, 0, 0, 0, 0, 0;
         0, 0, 0, 0, -1, 0;
         0, 0, 0, 0, 1, 0;
         0, -1, 0, 0, 0, 0;
         0, 1, 0, 0, 0, 0;
         0, 0, 0, 0, 0, -1;
         0, 0, 0, 0, 0, 1;
    ];

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
                [elNext, scoreOut, is_found] = steepestSlope(X, triples, frequencies, curEl, score, displacements, adopted, added, alpha, betta, gamma);
                if is_found == 0
                    break;
                end

                % traverse to this direction until minimum is reached
                while is_found
                    traversal_done = true;
                    elReserve = elNext;
                    [elNext, scoreOut, is_found] = steepestSlope(X, triples, frequencies, elNext, scoreOut, displacements, adopted, added, alpha, betta, gamma);
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








































% displacements = [
%     -1, 0, 0, 0, 0, 0;
%     -2, 0, 0, 0, 0, 0;
%      1, 0, 0, 0, 0, 0;
%      2, 0, 0, 0, 0, 0;
%      0, 0, 0, 0, -1, 0;
%      0, 0, 0, 0, -2, 0;
%      0, 0, 0, 0, 1, 0;
%      0, 0, 0, 0, 2, 0;
%      -1, 0, 0, 0, -1, 0;
%      1, 0, 0, 0, 1, 0;
%     ];
% 
% for i = 1:len % for each element
%     cur_el = adopted(i, :);
%     added = [];
%     
%     % left first
%     score1 = [0, 0];
%     cur = 0;
%     for j = 1:8
%         test = cur_el + displacements(j,:);
%         cur = cur + 1;
%         if testElement(test, 9)
%             score1(cur) = computeScore(X, frequencies, adopted, alpha, betta);
%         else
%             score1(cur) = score + 1;
%         end
%         
%         if cur == 2
%             % compare score1 elements to score
%             if max(score1) < score % we have to add the smallest element
%                 
%             end
%             score1 = [0, 0];
%             cur = 0;
%         end
%     end

        




