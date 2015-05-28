function selm = gpols_selection(popu,gap,pc,pm,tsels);
%Selection operator, mix the new generation
% selm = gpols_selection(popu,ggap,pc,pm,tsles)
%   selm <- matrix (n x 3)
%   popu -> population
%   ggap -> generation gap (0-1)
%   pc -> crossover probability (0-1)
%   pm -> mutation probability (0-1)
%   tsels -> selection (integer)
%
% Remark:
%   if tsels = 0 -> roulette wheel selection
%   if tsels = 1 -> total random selection
%   if tsels >= 2 -> tournament selection, tsels = tournament size
%  Example values: ggap = 0.8, pc = 0.7, pm = 0.3, tsels = 2
%  Columns of selm:
%   1: index of individual (first parent)
%   2: index of second ind. (second parant) (if crossover)
%   3: 0: direct rep., 1: crossover, 2: mutation
%

% (c) Janos Madar, University of Veszprem, 2005

%Begin
popun = popu.size;
selm = zeros(popun,3);

%Fitness values and sort
fit = zeros(1,popun);
for i = 1:popun,
  fit(i) = popu.chrom{i}.fitness;
end
if ~isempty(find(fit<0)),
  error('Every fitness must be > 0');
end
[fitsort,sortix] = sort(-fit);
fitsort = -fitsort./sum(-fitsort);
fitsum = cumsum(fitsort);
fitsum(end) = 1;  %avoid ~1E-16 error from representation

%Copy elite indv.s
i = floor((1-gap)*popun);
if i>=1 & i<=popun,
  selm(1:i,1) = sortix(1:i)';
  i = i+1;
else
  i = 1;
end
nn = i-1;

%New individuals
while nn<popun,
  if tsels > 0,
    j1 = tournament(fit,tsels);
  else
    j1 = roulette(fitsum,sortix);
  end
  %Select a method
  r = rand;
  if r<pc,
    %Crossover
    if tsels > 0,
      j2 = tournament(fit,tsels);
    else
      j2 = roulette(fitsum,sortix);
    end
    selm(i,1) = j1;
    selm(i,2) = j2;
    selm(i,3) = 1;
    i = i+1;
    nn = nn+2;
  elseif r<pc+pm,
    %Mutation
    selm(i,1) = j1;
    selm(i,3) = 2;
    i = i+1;
    nn = nn+1;
  else
    %Direct copy
    selm(i,1) = j1;
    selm(i,3) = 0;
    i = i+1;
    nn = nn+1;
  end
end
selm = selm(1:i-1,:);

%--------------------------------------------------------------
function j = tournament(fit,tsels);

n = length(fit);
jj = floor(rand(tsels,1)*n)+1;
[fitmax,maxix] = max(fit(jj));
j = jj(maxix);

%--------------------------------------------------------------
function j = roulette(fitsum,sortix);

v = find(fitsum >= rand(1,1));
j = sortix(v(1));
