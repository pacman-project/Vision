function tree = tree_getsubtree(treein,rix)
%Gets a subtree of a tree
% tree = tree_getsubtree(treein,rix)
%   tree <- the subtree
%   treein -> the source tree
%   rix -> subtree-index in source tree
%   symbols -> cell arrays of operator and terminator node strings
%

% (c) Janos Madar, University of Veszprem, 2005

%Begin
tree = treein;
if rix>treein.maxsize, return; end
%Get the subtree
vin  = [rix];
vout = [1];
iix = 1;
while vin(iix)<=treein.maxsize,
  tree.nodetyp(vout(iix)) = treein.nodetyp(vin(iix));
  tree.node(vout(iix)) = treein.node(vin(iix));
  vin = [vin vin(iix)*2 vin(iix)*2+1];
  vout = [vout vout(iix)*2 vout(iix)*2+1];
  iix = iix+1;
end
