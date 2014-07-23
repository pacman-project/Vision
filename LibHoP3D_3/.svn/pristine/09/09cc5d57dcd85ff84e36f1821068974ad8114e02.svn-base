% this is to compute curvedness of the element


function [ curvedness ] = compute3Curvedness(el)

    central = el(3:4);
    left = el(1:2);
    right = el(5:6);
    
    % distance between left and right is to prevent wired elements
    
    % first term is curvedness, second term is skewness
    
    % curvednessD = sqrt(sum((central - left).^2)) + sqrt(sum((central - right).^2));
    curvednessD = 0;
    
    skewnessD = sqrt( (left(2) - right(2))^2 + (central(2) - right(2))^2  + (central(2) - left(2))^2);
    curvedness = skewnessD + curvednessD;
end

