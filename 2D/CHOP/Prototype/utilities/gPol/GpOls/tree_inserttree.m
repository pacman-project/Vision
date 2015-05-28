function tree = tree_inserttree(treesrc,treedst,dstix,nvar)
%Inserts a tree into a tree (as a subtree)
% tree = tree_inserttree(treesrc,treedst,dstix,nvar)
%   tree <- result
%   treesrc -> source tree, it will be inserted into destination tree
%   treedst -> destination tree
%   dstix -> insert-index in destination tree
%   nvar -> number of terminator node types (variables)
%

% (c) Janos Madar, University of Veszprem, 2005

%Begin
tree = treedst;
if dstix>treedst.maxsize, return; end
%Insert
vin = [1];
vout = [dstix];
iix = 1;
while vin(iix)<=treesrc.maxsize & vout(iix)<=tree.maxsize,
  tree.nodetyp(vout(iix)) = treesrc.nodetyp(vin(iix));
  tree.node(vout(iix)) = treesrc.node(vin(iix));
  vin = [vin vin(iix)*2 vin(iix)*2+1];
  vout = [vout vout(iix)*2 vout(iix)*2+1];
  iix = iix+1;
end
%Repeair the tail
for i = (tree.maxsize+1)/2:tree.maxsize,
  if tree.nodetyp(i)==1,
    tree.nodetyp(i) = 2;
    tree.node(i) = floor(nvar*rand)+1;
  end
end
