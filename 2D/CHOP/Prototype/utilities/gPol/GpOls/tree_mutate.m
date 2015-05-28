function tree = tree_mutate(treein,symbols);
%Mutates a tree (mutates one randomly selected node)
% tree = tree_mutate(treein,symbols)
%   tree <- the output tree
%   treein -> the input tree
%   symbols -> cell arrays of operator and terminator node strings
%

% (c) Janos Madar, University of Veszprem, 2005

%Begin
tree = treein;
nn = [length(symbols{1}), length(symbols{2})];
%Mutate one node
[n,v] = tree_size(tree);
i = v(floor(rand*(length(v))+1));
if i<(tree.maxsize+1)/2 & rand<0.5,
  [tree.nodetyp(i) tree.node(i)] = tree_genrndsymb((tree.nodetyp(i)==1),nn);
else
  while tree.node(i)==treein.node(i) & tree.nodetyp(i)==treein.nodetyp(i),
    [tree.nodetyp(i) tree.node(i)] = tree_genrndsymb((tree.nodetyp(i)~=1),nn);
  end
end

%------------------------------------------------------------------
function [nodetyp,node] = tree_genrndsymb(p0,nn)
%Generate a random symbol (tarminate or operate)
%  [nodetyp,node] = tree_genrndsymb(p0,nn)
%   nodetyp,node <- results
%   p0 -> probability of terminate node
%   nn -> vector [number of operators, variables]
%

if rand<p0,
  nodetyp = 2;
else
  nodetyp = 1;
end
node = floor(nn(nodetyp)*rand)+1;