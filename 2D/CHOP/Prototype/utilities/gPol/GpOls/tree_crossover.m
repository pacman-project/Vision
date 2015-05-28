function [tree1,tree2] = recombinate_trees(treein1,treein2,mode,symbols);
%Recombinates two trees
% [tree1,tree2] = recombinate_trees(treein1,treein2,mode,symbols)
%   tree1,tree2 <- two childs
%   treein1,treein2 -> two parents
%   mode -> 1: one-point- 2: two-point crossover
%   symbols -> cell arrays of operator and terminator node strings
%

% (c) Janos Madar, University of Veszprem, 2005

%Begin
tree1 = treein1;
tree2 = treein2;
nn = [length(symbols{1}), length(symbols{2})];
%Calculate indexes
switch mode,
  case 1,
    [n,v1] = tree_size(tree1);
    [n,v2] = tree_size(tree2);
    n = max([tree1.maxsize, tree2.maxsize]);
    dummy1 = zeros(n,1);
    dummy2 = zeros(n,1);
    dummy1(v1) = 1;
    dummy2(v2) = 1;
    v = find((dummy1+dummy2)==2);
    ix1 = v(floor(rand*(length(v))+1));
    ix2 = ix1;
  case 2,
    [n,v] = tree_size(tree1);
    ix1 = v(floor(rand*(length(v))+1));
    [n,v] = tree_size(tree2);
    ix2 = v(floor(rand*(length(v))+1));
  otherwise,
    return;
end
%Repleace subtrees (recombinate)
sub1 = tree_getsubtree(treein1,ix1);
sub2 = tree_getsubtree(treein2,ix2);
tree1 = tree_inserttree(sub2,tree1,ix1,nn(2));
tree2 = tree_inserttree(sub1,tree2,ix2,nn(2));
