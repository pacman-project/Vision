function popu = gpn_evaluate(popuin,ixs,X,Y,Q,optv);
%Evaluates individulas and identificates their linear parameters
% popu = gpols_evaluate(popuin,ixs,X,Y,Q,optv)
%   popu <- result (population)
%   popuin -> input population
%   ixs -> vector of indexes of individulas to evaluate
%   X -> Regression matrix without bias (!)
%   Y -> Output vector
%   Q -> Weighting vector (set empty if not used)
%   optv -> evaluation parameters (set empty for default)
%    [optv(1) optv(2)]: a1, a2 tree-size penalty parameters (default: 0,0)
%    optv(3): OLS treshold value, range: 0-1 (default: 0)
%    optv(4): if == 1 then polynomial evaluation else normal (default: 0)
%
%Remark:
%  The loss function is:
%    mse = (1/n)*E'*diag(Q)*E
%   where n is the number of data points, E is the error vector (Y-Ym),
%   diag(Q) is the weighting diagonal matrix (e.g. inverse of noise
%   covariance matrix or relative weights for relative mse)
%  The Y is the given output vector and the Em is the estimated output:
%    Y = Ym + E, Ym = [X ones(size(X,1),1)]*Theta
%   where X is the regression matrix, Theta is the linear parameter vector
%   and Theta(end) = bias
%
%  The ixs contains the indexes of individuals must be evaluated.
%   If you want evaluate every individudal then set ixs = [1:popu.size]
%
%  The applied tree-size penalty function:
%    penalty = 1/(1+exp(a1*(L-a2)))
%   where a1 and a2 the penalty function parameters, and L is the size of
%   the tree (number of nodes). If a1=0 or a2=0 then penalty = 1
%  Then the fitness is calucalted as
%    fitness = penalty * corrcoef(sqrt(diag(Q))*Y,sqrt(diag(Q))*Ym)
%
%  The OLS treshold value determines the minimum allowed error reduction
%  ratio (err). If a sub-tree has smaller err then it will be eliminated.
%
%  The polynomial evaluation means that every non-polynomial tree will be
%  transformed to polynomial before evaluation.
%

% (c) Janos Madar, University of Veszprem, 2005

%Output
popu = popuin;

%Options and parameters
if isempty(optv),
  optv = zeros(1,4);
end
a1 = optv(1);
a2 = optv(2);
olslimit = optv(3);
polye = optv(4);

%WLS matrices
if isempty(Q),
  Q = ones(size(Y,1),1);
end
Q = diag(Q);
X = sqrt(Q)*X;
Y = sqrt(Q)*Y;

%Symbolum list
for i = 1:length(popu.symbols{1}),
  s = popu.symbols{1}{i}(1);
  if s=='*' | s=='/' | s=='^' | s=='\',
    symbols{1}{i} = strcat('.',popu.symbols{1}{i});
  else
    symbols{1}{i} = popu.symbols{1}{i};
  end
end
for i = 1:size(X,2),
  symbols{2}{i} = sprintf('X(:,%i)',i);
end

%MAIN loop
for j = ixs,

  %Get the tree
  tree = popu.chrom{j}.tree;

  %Exhange '+' -> '*' under non-'+' (polynom-operation)
  if (polye == 1),
    tree = polytree(tree);
  end

  %Collect the '+ parts'
  [vv,fs,vvdel] = fsgen(tree,symbols);

  %Prune redundant parts
  tree = prunetree(vvdel,tree,symbols);

  %Collect the '+ parts'
  [vv,fs,vvdel] = fsgen(tree,symbols);
  if ~isempty(vvdel),
    error('Fatal error: redundant strings after deleting');
  end

  %OLSQ
  [vfsdel,err] = gpols_olsq(fs,X,Y,olslimit);
  tree.err = err;
  vvdel = vv(vfsdel);

  %Prune redundant parts
  tree = prunetree(vvdel,tree,symbols);

  %Collect the '+ parts'
  [vv,fs,vvdel] = fsgen(tree,symbols);
  if ~isempty(vvdel),
    error('Fatal error: redundant strings after deleting');
  end

  %LSQ
  [mse,cfsq,theta] = gpols_lsq(fs,X,Y);
  fit = cfsq;

  %Tree-size penalty
  if a1~=0 & a2~=0,
    Sl = tree_size(tree);
    fit = fit / (1+exp(a1*(Sl-a2)));
  end

  %Chrom
  popu.chrom{j}.tree = tree; %write back the tree
  popu.chrom{j}.mse = mse;
  popu.chrom{j}.fitness = fit;
  popu.chrom{j}.tree.param(1:length(theta)) = theta;
  popu.chrom{j}.tree.paramn = length(theta);

end

%---------------------------------------------------------------
function [tree] = polytree(treein);
tree = treein;
v = [1];
vv = [];
i = 1;
while i <= length(v),
  ii = v(i);
  if tree.nodetyp(ii)==1 & tree.node(ii)==1,
    v = [v, ii*2, ii*2+1];
  else
    vv = [vv, ii];
  end
  i = i+1;
end
for ii = [vv],
  v = [ii];
  i = 1;
  while i <= length(v),
    if tree.nodetyp(v(i))==1,
      if tree.node(v(i))==1,
        tree.node(v(i)) = 2;
      end
      if v(i)*2+1 <= tree.maxsize,
        v = [v, v(i)*2, v(i)*2+1];
      end
    end
    i = i+1;
  end
end

%---------------------------------------------------------------
function [vv,fs,vvdel] = fsgen(tree,symbols);
%Search the '+ parts'
v = [1];
vv = [];
i = 1;
while i <= length(v),
  ii = v(i);
  if tree.nodetyp(ii)==1 & tree.node(ii)==1,
    v = [v, ii*2, ii*2+1];
  else
    vv = [vv, ii];
  end
  i = i+1;
end
fs = [];
i = 1;
for ii = [vv],
  fs{i} = strcat('(',tree_stringrc(tree,ii,symbols),')');
  i = i+1;
end
%Search the redundant '+ parts'
vvdel = [];
vvv = [];
i = 1;
while i <= length(fs),
  ok = 0;
  ii = 1;
  while ii<i & ok==0,
    ok = strcmp(fs{i},fs{ii});
    ii = ii+1;
  end
  if ok==1,
    vvdel = [vvdel, vv(i)];
  else
    vvv = [vvv, i];
  end
  i = i+1;
end

%---------------------------------------------------------------
function tree = prunetree(vvdel,treein,symbols);
%Delete subtrees
nn = [length(symbols{1}), length(symbols{2})];
tree = treein;
n = floor(tree.maxsize/2);
tree.nodetyp(vvdel) = 0;
ok = 1;
while ok,
  ok = 0;
  i = 1;
  while i<=n & ok==0,
    if (tree.nodetyp(i)==1) & (tree.nodetyp(i*2)==0 | tree.nodetyp(i*2+1)==0),
      ok = 1;
      if tree.nodetyp(i*2)==0 & tree.nodetyp(i*2+1)==0,
        tree.nodetyp(i*2) = treein.nodetyp(i*2);
        tree.nodetyp(i*2+1) = treein.nodetyp(i*2+1);
        tree.nodetyp(i) = 0;
      elseif tree.nodetyp(i*2)==0,
        tree.nodetyp(i*2) = treein.nodetyp(i*2);
        subtree = tree_getsubtree(tree,i*2+1);
        tree = tree_inserttree(subtree,tree,i,nn(2));
      else
        tree.nodetyp(i*2+1) = treein.nodetyp(i*2+1);
        subtree = tree_getsubtree(tree,i*2);
        tree = tree_inserttree(subtree,tree,i,nn(2));
      end
    elseif (tree.nodetyp(i*2)==0 | tree.nodetyp(i*2+1)==0),
      ok = 1;
      if tree.nodetyp(i*2)==0,
        tree.nodetyp(i*2) = treein.nodetyp(i*2);
      end
      if tree.nodetyp(i*2+1)==0,
        tree.nodetyp(i*2+1) = treein.nodetyp(i*2+1);
      end
    end
    i = i+1;
  end
end
