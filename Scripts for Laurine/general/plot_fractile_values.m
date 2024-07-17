function [ horizontal_lines, vertical_lines ] = plot_fractile_values( xmin, xmax, ymin, P_fractiles, x_fractiles, show1, show2 )
% plotting the fractile values


horizontal_lines =[];
vertical_lines =[];

for i = 1:length(P_fractiles)

horizontal_lines(i) = line('XData', [xmin xmax], 'YData', [P_fractiles(i) P_fractiles(i)], 'Color',[0 0 0]+0.05*10,'HandleVisibility','off');      
vertical_lines(i) = line('XData', [x_fractiles(i) x_fractiles(i)], 'YData', [ymin P_fractiles(i)], 'LineWidth', 0.01, 'Color', [0 0 0]+0.05*10, 'LineStyle','--', show1, show2);

end



end

