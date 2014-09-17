%> Name: fastintersect
%>
%> Description: Provides a fast intersection function between two sets of
%> positive integers. (Code found on the internet)
%>
%> @param A First set to intersect.
%> @param B Second set to intersect.
%>
%> @retval C Intersection set.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 17.02.2014
function C = fastintersect(A,B)
    if ~isempty(A)&&~isempty(B)
       P = zeros(1, max(max(A),max(B)) ) ;
       P(A) = 1;
       C = B(logical(P(B)));
    else
        C = [];
    end
end