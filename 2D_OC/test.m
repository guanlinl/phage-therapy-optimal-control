%%%%%%%%%%%%%%%%%%%%%
% Put legend in its own subplot
%%%%%%%%%%%%%%%%%%%%%
%
%   NOTE: The legend will not resize so set the size and then run it all
   h = figure(1); clf
%   
%   % Make the main plot in 1
   ax1 = subplot(1,2,1);
   plot(rand(10,3))
   lgh = legend('a','b','c');	% Set your legend
%   
%   % Put the legend in ax2
   ax2 = subplot(1,2,2);
   set(lgh,'position',get(ax2,'position'));
   axis(ax2,'off');
   title('legend')
%
%
