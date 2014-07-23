% function selects steepest slope

function [el_out, scoreOut, is_found] = steepestSlope4(X, triples4, frequencies, curEl, score, triples3Out, n3Clusters, adopted, added, ...
                                                        alpha, betta, gamma, maxDist, distThresh, distXtoInitAdopted)

%     str = ['Looking for steepest slope ', num2str(curEl)];
%     disp(str);

    % parts of the curEl
    el_top = curEl(1:6);
    el_middle = curEl(7:12);
    el_bottom = curEl(13:18);
    
    points = [el_top; el_middle; el_bottom];
    
    % compute distances between all elements of the third layer
    distances3 = Isodata_distances(points, triples3Out, 3, n3Clusters, false, false);
    
    % find at least numDisp closest elements
    a1 = sort(distances3(:), 'ascend');
    a1 = a1(4:end); % first three are zeros
    
    r = [];
    c = [];
    ind = 1;
    while a1(ind) < maxDist %(length(r) < numDisp)
        [r1, c1] = find(distances3 == a1(ind));
        r = [r;r1];
        c = [c;c1];
        ind = ind + 1;        
    end
    
    lenDisp = length(r);
    displacements = zeros(lenDisp, 18);
    
    for i = 1: lenDisp;
        start = (r(i)-1)*6 + 1;
        el = triples3Out(c(i),:);
        displacements(i,:) = curEl;
        displacements(i, start:start+5) = el;
    end
    
    is_found = false;
    el_out = curEl;
    scoreOut = score;
    
    for j = 1:lenDisp
        el = displacements(j,:);
        
        % then we test whether the element's frequency is larger than zero
        [topIndex, ~] = compute3elementIndex(el(1:6), triples3Out, n3Clusters);
        [middleIndex, ~] = compute3elementIndex(el(7:12), triples3Out, n3Clusters);
        [bottomIndex, ~] = compute3elementIndex(el(13:18), triples3Out, n3Clusters);
        
        if triples4(topIndex, middleIndex, bottomIndex) > 0 % if this element appeares at least ones in the data
            scoreTrial = computeScore(X, frequencies, adopted, [added; el], alpha, betta, distThresh, distXtoInitAdopted);
            scoreTrial = scoreTrial + gamma * compute4Curvedness(el);
            if scoreTrial < scoreOut % there is a slope
                scoreOut = scoreTrial;
                el_out = el;
                is_found = true;
            end
        end
    end

end

