function bubble_contours = get_annotation_contour(mst_tree, node_positions, node_size, bubbles, r, resolution)

figure(100); draw_SPADE_tree(mst_tree, node_positions, node_size, [], [], 0, 0);
canvas_axis = axis; close(100);

dist = zeros(length(canvas_axis(1):(canvas_axis(2)-canvas_axis(1))/resolution:canvas_axis(2))*length(canvas_axis(3):(canvas_axis(4)-canvas_axis(3))/resolution:canvas_axis(4)),size(node_positions,2));
grid_position = zeros(size(dist,1),2);
% resolution = 100;
delta_x = (canvas_axis(2)-canvas_axis(1))/resolution;
delta_y = (canvas_axis(4)-canvas_axis(3))/resolution;
counter = 1;
for x=canvas_axis(1):delta_x:canvas_axis(2)
    for y=canvas_axis(3):delta_y:canvas_axis(4)
        dist(counter,:) = sqrt((node_positions(1,:)-x).^2 + (node_positions(2,:)-y).^2);
        grid_position(counter,:) = [x,y];
        counter = counter + 1;
    end
end

% figure(100); draw_SPADE_tree(mst_tree, node_positions, node_size, [], [], 0, 0);
% figure(100); hold on;
for i=1:length(bubbles)
    indicators = (sum(dist(:,bubbles{i})<r(i),2)~=0) & (min(dist(:,bubbles{i})')'+r(i)/10<min(dist(:,setdiff(1:end,bubbles{i}))')');
    outer_points = get_outer_points(grid_position, indicators, delta_x, delta_y);
%     line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)])
%     for i=1:size(outer_points,1), plot(outer_points(i,1),outer_points(i,2),'g.'); pause(0.05); end
%     for i=1:size(outer_points,1)-1, line(outer_points(i:i+1,1),outer_points(i:i+1,2)); pause(0.05); end
    bubble_contours{i} = outer_points;
end

return




function outer_points = get_outer_points(grid_position, indicators, delta_x, delta_y)

patch_squareform = reshape(indicators,sqrt(length(indicators)), sqrt(length(indicators)));
contour_squareform = zeros(size(patch_squareform));

for i=1:size(patch_squareform,1)
    ind_tmp = (patch_squareform(i,:)~=[-1,patch_squareform(i,1:end-1)] | patch_squareform(i,:)~=[patch_squareform(i,2:end),-1]) & patch_squareform(i,:)==1;
    contour_squareform(i,ind_tmp)=1;
end
for j=1:size(patch_squareform,2)
    ind_tmp = (patch_squareform(:,j)~=[-1;patch_squareform(1:end-1,j)] | patch_squareform(:,j)~=[patch_squareform(2:end,j);-1]) & patch_squareform(:,j)==1;
    contour_squareform(ind_tmp,j)=1;
end

contour_squareform = contour_squareform(:);
outer_points = grid_position(contour_squareform==1,:);

% order them
ordered_points = outer_points(1,:); 
flag_outer_points = zeros(size(outer_points,1),1);flag_outer_points(1)=1;
while sum(flag_outer_points)~=length(flag_outer_points)
    dist = sum((repmat(ordered_points(end,:),size(outer_points,1),1) - outer_points).^2,2);
    dist(flag_outer_points==1) = max(dist) + 1;
    [dummy,I] = min(dist);
    if dummy > sqrt(delta_x^2+delta_y^2)*2;
        flag_outer_points(I)=1; continue;
    end
    ordered_points = [ordered_points; outer_points(I,:)]; flag_outer_points(I)=1;
end
outer_points = ordered_points;