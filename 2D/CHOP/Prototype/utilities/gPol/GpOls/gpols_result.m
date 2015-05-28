function [sout,tree] = gpols_result(popu,info);
%Gets information string about the best solution of a population
%  [sout,tree] = gpols_result(popu,info);
%   sout <- text (string)
%   tree <- the best solution
%   popu -> population structure
%   info -> info mode (1,2)
%

% (c) Janos Madar, University of Veszprem, 2005

if info == 0,
  sout = sprintf('Iter \t Fitness \t Solution');
  return;
end

best = popu.chrom{1}.fitness;
bestix = 1;
for i = 1:popu.size,
  if popu.chrom{i}.fitness > best,
    best = popu.chrom{i}.fitness;
    bestix = i;
  end
end
tree = popu.chrom{bestix}.tree;

if info == 1,
  sout = sprintf('%3i. \t %f',popu.generation,best);
  s = tree_stringrc(tree,1,popu.symbols);
  sout = sprintf('%s \t %s',sout,s);
  return;
end

if info == 2,
  sout = sprintf('fitness: %f,  mse: %f',best,popu.chrom{bestix}.mse);
  [vv,fs] = fsgen(tree,popu.symbols);
  for i = 1:length(fs),
    sout = sprintf('%s\n %f * %s +',sout,tree.param(i),fs{i});
  end
  sout = sprintf('%s\n %f',sout,tree.param(i+1));
  return;
end

sout = '???';

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

  