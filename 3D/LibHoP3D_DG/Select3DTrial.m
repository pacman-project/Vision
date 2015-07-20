 h = surf(peaks);
 zoom(10);
 disp('Click anywhere on the surface, then hit return')
 pause
 [p v vi face facei] = select3d;
 marker1 = line('xdata',p(1),'ydata',p(2),'zdata',p(3),'marker','o',...
                'erasemode','xor','markerfacecolor','k');
 marker2 = line('xdata',v(1),'ydata',v(2),'zdata',v(3),'marker','o',...
                'erasemode','xor','markerfacecolor','k');
 marker2 = line('erasemode','xor','xdata',face(1,:),'ydata',face(2,:),...
                'zdata',face(3,:),'linewidth',10);
 disp(sprintf('\nYou clicked at\nX: %.2f\nY: %.2f\nZ: %.2f',p(1),p(2),p(3)'))
 disp(sprintf('\nThe nearest vertex is\nX: %.2f\nY: %.2f\nZ: %.2f',v(1),v(2),v(3)'))

