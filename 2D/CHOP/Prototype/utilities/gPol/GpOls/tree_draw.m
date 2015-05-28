function tree_draw(tree,symbols);
%Creates a figure and draws the tree
%  tree_draw(tree,symbols)
%   tree -> the tree
%   symbols -> cell arrays of operator and terminator node strings
%

% (c) Janos Madar, University of Veszprem, 2005


figure1 = figure('Color',[0.8 0.8 0.8]);
set(figure1,'DefaultAxesFontName','times');
set(figure1,'DefaultTextFontName','times');
set(figure1,'DefaultAxesFontSize',8);
set(figure1,'DefaultTextFontSize',8);
ax1 = axes('Position',[0 0 1 1]);
subplot(ax1);
axis([0 1 0 1]);
axis off
 
[n,v] = tree_size(tree);
iv = 1;
maxrow = floor(log(max(v))/log(2)+1);
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
 %Put text
 j = v(iv);
 s = symbols{tree.nodetyp(j)}{tree.node(j)};
 set(label,'String',s);
 %Next
 iv = iv+1;
end;

