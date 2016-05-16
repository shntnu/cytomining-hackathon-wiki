function [highest_density_position] = find_highest_density_position(cell_positions, local_density)

[F,XI,hx] = ksdensity(cell_positions(1,:)); 
[F,YI,hy] = ksdensity(cell_positions(2,:)); 
grid_points = [repmat(XI(:)',1,length(YI));reshape(repmat(YI(:)',length(XI),1),1,length(XI)*length(YI))];
window_x = abs(XI(2)-XI(1))*1.1;
window_y = abs(YI(2)-YI(1))*1.1;
freq = zeros(2,length(XI)*length(YI))+NaN;
for i=1:size(grid_points,2)
    freq(1,i) = sum(exp(-(cell_positions(1,:) - grid_points(1,i)).^2/(2*hx^2)).*exp(-(cell_positions(2,:) - grid_points(2,i)).^2/(2*hy^2))/(sqrt(2*pi*hx^2*2*pi*hy)));
    freq(2,i) = sum((exp(-(cell_positions(1,:) - grid_points(1,i)).^2/(2*hx^2)).*exp(-(cell_positions(2,:) - grid_points(2,i)).^2/(2*hy^2))/(sqrt(2*pi*hx^2*2*pi*hy))).*local_density);
end
freq(1,:) = freq(1,:)/size(cell_positions,2);
freq(2,:) = freq(2,:)/sum(local_density);

[~,ind] = max(freq(2,:));
% contour plot of original data
highest_density_position = grid_points(:,ind);

