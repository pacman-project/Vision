%> Name: calculateCoveredPixels
%>
%> Description: Given a set of leaf nodes and program options, this
%> function calculates whch pixels of every image is covered by which leaf
%> node. The info is kept inside a cell array, and returned in 'coveredNodes'
%> variable. Each cell includes the linear indices of covered pixels. 
%>
%> @param leafNodes Leaf node array of the form 
%> [labelId, posX, posY, imageId; ...].
%> @param options program options.
%>
%> @retval coveredNodes The cell array that includes a set of covered
%> pixels' linear indices for every node.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 16.09.2015
function [ output_args ] = calculateCoveredPixels( leafNodes, options )
%CALCULATECOVEREDPIXELS Summary of this function goes here
%   Detailed explanation goes here


end

