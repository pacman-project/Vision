%> Name: fastsortedunique
%>
%> Description: As it is clear from the name, calculates the unique set of a
%> given sorted arrayn of elements.
%>
%> @param vector The element array.
%>
%> @retval vector Unique set of elements.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.09.2014
function [ vector ] = fastsortedunique( vector )
   vector = vector([true;diff(vector(:))>0]);
end

