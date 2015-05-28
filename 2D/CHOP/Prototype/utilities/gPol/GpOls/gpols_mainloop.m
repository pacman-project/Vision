function [popu,evnum] = gpols_mainloop(popuin,X,Y,Q,opt)
%Run one evolutionary loop, makes the next generation
% [popu,evnum] = gpols_mainloop(popuin,X,Y,Q,opt)
%   popu <- next generation of the population
%   evnum <- number of fun.evaulation (usually number of new individuals)
%   popuin -> the population
%   opt -> options vector, GPOLS-parameters
%   X,Y,Q -> input, output and weighting matrices (see gpols_evaluate)
%
% Remark:
%     opt(1): ggap, generation gap (0-1)
%     opt(2): pc, probability of crossover (0-1)
%     opt(3): pm, probability of mutation (0-1)
%     opt(4): selection type (integer, see gpols_selection)
%     opt(5): rmode, mode of tree-recombination (1 or 2)
%     opt(6): a1, first penalty parameter
%     opt(7): a2, second penalty parameter (0 if there is not penalty)
%     opt(8): OLS treshhold real 0-1 or integer >= 2
%     opt(9): if == 1 -> polynomial evaluation
%     opt(10): if == 1 -> evaluate all indv.s not only new offsprings
%

% (c) Janos Madar, University of Veszprem, 2005

popun = popuin.size;
ggap = opt(1);
pc = opt(2);
pm = opt(3);
tsels = opt(4);
rmode = opt(5);

%Selection
selm = gpols_selection(popuin,ggap,pc,pm,tsels);

%New generation
popu = popuin;
newix = [];
nn = 1;
for i=1:size(selm,1),
  m = selm(i,3);
  %*** Crossover ***
  if m==1,
    p1 = selm(i,1);
    p2 = selm(i,2);
    popu.chrom{nn} = popuin.chrom{p1};
    popu.chrom{nn}.fitness = -1;
    if nn+1<popu.size, 
      popu.chrom{nn+1} = popuin.chrom{p2};
      popu.chrom{nn+1}.fitness = -1;
    end
    %recombinate trees
    tree1 = popuin.chrom{p1}.tree;
    tree2 = popuin.chrom{p2}.tree;
    [tree1,tree2] = tree_crossover(tree1,tree2,rmode,popu.symbols);
    popu2.chrom{nn}.tree = tree1;
    if nn+1<=popu.size,
      popu.chrom{nn+1}.tree = tree2;
    end
    %remember the new individulas
    newix = [newix nn];
    nn = nn+1;
    if nn<=popu.size,
      newix = [newix nn];
      nn = nn+1;
    end
  %*** Mutation ***
  elseif m==2,
    p1 = selm(i,1);
    popu.chrom{nn} = popuin.chrom{p1};
    popu.chrom{nn}.fitness = -1;
    %muatate tree
    tree1 = popu.chrom{p1}.tree;
    tree1 = tree_mutate(tree1,popu.symbols);
    popu.chrom{nn}.tree = tree1;
    %remember the new individual
    newix = [newix nn];
    nn = nn+1;
  %*** Direct copy *** 
  else
    p1 = selm(i,1);
    popu.chrom{nn} = popu.chrom{p1};
    nn = nn+1;
  end
end

%if opt(10)==1 -> evaluate all indv.s
if length(opt)>9 & opt(10)==1,
  newix = [1:popu.size];
end

%Evaulate new indv.s
popu = gpols_evaluate(popu,newix,X,Y,Q,opt(6:9));
popu.generation = popu.generation+1;
evnum = length(newix);
