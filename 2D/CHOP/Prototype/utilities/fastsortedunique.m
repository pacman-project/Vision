function [ vector ] = fastsortedunique( vector )
   vector = vector([true;diff(vector(:))>0]);
end

