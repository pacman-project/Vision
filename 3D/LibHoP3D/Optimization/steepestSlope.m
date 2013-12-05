% function selects steepest slope

function [el_out, scoreOut, is_found] = steepestSlope(X, triples, frequencies, curEl, score, displacements, adopted, added, alpha, betta, gamma)

%     str = ['Looking for steepest slope ', num2str(curEl)];
%     disp(str);

    [rD, ~] = size(displacements);
    is_found = false;
    el_out = curEl;
    scoreOut = score;
    
    for j = 1:rD
        displ = displacements(j,:);
        el = curEl + displ;
        [el2Ind, is_ok] = test3Element(el, 9);
        % then we test whether the element's frequency is larger than zero
        
        if is_ok && triples(el2Ind(1), el2Ind(2), el2Ind(3))>0
            scoreTrial = computeScore(X, frequencies, [adopted; added; el], alpha, betta);
            scoreTrial = scoreTrial + gamma * compute3Curvedness(el);
            if scoreTrial < scoreOut % there is a slope
                scoreOut = scoreTrial;
                el_out = el;
                is_found = true;
            end
        end
    end

end

