function tree = tree_genuser(treein,symbols);
%Generates a tree by user
%  tree = genuser_tree(treein,symbols);
%   tree <- the result
%   treein -> begining tree
%   symbols -> cell arrays of operator and terminator node strings
%
%  Remark: You must initialize the tree, before call this function
%

% (c) Janos Madar, University of Veszprem, 2005

tree = treein;

figure1 = figure('Color',[0.8 0.8 0.8]);
set(figure1,'DefaultAxesFontName','times');
set(figure1,'DefaultTextFontName','times');
set(figure1,'DefaultAxesFontSize',8);
set(figure1,'DefaultTextFontSize',8);
ax1 = axes('Position',[0 0 1 1]);
subplot(ax1);
axis([0 1 0 1]);
axis off
v = [1];
iv = 1;
maxrow = floor(log(tree.maxsize)/log(2)+1);
dy = 0.95/maxrow;
while iv<=length(v),
  %Draw
  j = v(iv);
  row = floor(log(j)/log(2)+1);
  col = j - 2^(row-1);
  x1 = (2*col+1)/(2^row);
  y1 = 1-row*dy;
  label = text(x1,y1,'(?)');
  if row>1,
    j = floor(j/2);
    row2 = floor(log(j)/log(2)+1);
    col2 = j - 2^(row2-1);
    x2 = (2*col2+1)/(2^row2);
    y2 = 1-row2*dy;
    line([x1;x2],[y1+dy/8;y2-dy/8]);
  end
  %Get symbol
  j = v(iv);
  ok = 1;
  while ok,
    s = input(':','s');
    i = 1;
    while i<=length(symbols{1}) & ~strcmp(symbols{1}{i},s), i=i+1; end
    if i<=length(symbols{1}) & j<(tree.maxsize+1)/2,
      v = [v, j*2, j*2+1];
      tree.nodetyp(j) = 1;
      tree.node(j) = i;
      break;
    end
    i = 1;
    while i<=length(symbols{2}) & ~strcmp(symbols{2}{i},s), i=i+1; end
    if i<=length(symbols{2}),
      tree.nodetyp(j) = 2;
      tree.node(j) = i;
      break;
    end
  end
  %Draw symbol
  set(label,'String',s);
  %Next
  iv = iv+1;
end

