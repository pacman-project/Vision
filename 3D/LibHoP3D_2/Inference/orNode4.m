% this is an imitation of OR-nodes used in the indexin process for or nodes
% el = [top, central, bottom] which are 3layer elements

function newEl = orNode4(el, table3, triples3Out)

    derivatives = zeros(3,6);
    
    adderTop = zeros(1,6);
    adderBottom = zeros(1,6);
    
    
    for i = 1:3
        curEl = el(i);
        % extract derivatives for this element
        derivatives(i,:) = triples3Out(curEl, :);
    end
    
    is_changed = false;    
    for i = 2:2:6 % y-s
        
        if derivatives(2,i) - derivatives(1,i) == 2 || derivatives(2,i) - derivatives(1,i) == 3
            adderTop(i) = derivatives(2,i) - derivatives(1,i) - 1;
        elseif derivatives(2,i) - derivatives(1,i) == -2 || derivatives(2,i) - derivatives(1,i) == -3
            adderTop(i) = derivatives(2,i) - derivatives(1,i) + 1;
        end
        
        if derivatives(2,i) - derivatives(3,i) == 2 || derivatives(2,i) - derivatives(3,i) == 3
            adderBottom(i) = derivatives(2,i) - derivatives(3,i) - 1;
        elseif derivatives(2,i) - derivatives(3,i) == -2 || derivatives(2,i) - derivatives(3,i) == -3
            adderBottom(i) = derivatives(2,i) - derivatives(3,i) + 1;
        end
        
    end
    
    if adderTop(2) == adderTop(4) && adderTop(4) == adderTop(6)
        derivatives(1,:) = derivatives(1,:) + adderTop;
        is_changed = true;
    end
    if adderBottom(2) == adderBottom(4) && adderBottom(4) == adderBottom(6)
        derivatives(3,:) = derivatives(3,:) + adderBottom;
        is_changed = true;
    end
    
    % recompute els
    newEl = el;
    
    if is_changed
        for i = 1:2:3  % 1 and 3
            cur = derivatives(i,:);
            if table3(cur(1),cur(2),cur(3),cur(4),cur(5),cur(6)) ~= 0 % this element exists
                newEl(i) = table3(cur(1),cur(2),cur(3),cur(4),cur(5),cur(6));
            end
        end
    end

end

