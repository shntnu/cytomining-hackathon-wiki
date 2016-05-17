function [tree_bubble_contours] = draw_SPADE_tree(mst_tree, node_positions, node_size, node_data, data, show_color, show_annotation,is_selected, color_scheme, tree_annotations,tree_bubble_contours)

if ~exist('tree_bubble_contours')
    tree_bubble_contours = cell(size(tree_annotations));
end

% clear the plot
hold off; plot(0,0,'visible','off'); hold on;
% draw edges
adj = mst_tree;
coeff = node_positions;
pairs = SPADE_find_matrix_big_element(triu(adj,1),1);
for k=1:size(pairs,1), line(coeff(1,pairs(k,:)), coeff(2,pairs(k,:)),'color','g'); end

% % show annotation contour
% show_annotation = 0;
switch show_annotation % 0 = no show; 1 = show all; 2 = show selected
    case 0 
        'no show skip bubble';
    case 1 

        % compute bubble contours
        bubbles = tree_annotations;
        annotation_bubble_size = 5;
        r = ones(1,length(bubbles))*annotation_bubble_size; 
        resolution = 100;

        % % control axis limits
        axis_lims = reshape([-max(abs(coeff)'); +max(abs(coeff)')], 1, 4)*1.1;
        for i=1:4
            if abs(axis_lims(i))<55
                if axis_lims(i)>=0
                    axis_lims(i)=55;
                else
                    axis_lims(i)=-55;
                end
            end
        end
        canvas_axis = axis_lims; 

        dist = zeros(length(canvas_axis(1):(canvas_axis(2)-canvas_axis(1))/resolution:canvas_axis(2))*length(canvas_axis(3):(canvas_axis(4)-canvas_axis(3))/resolution:canvas_axis(4)),size(node_positions,2));
        grid_position = zeros(size(dist,1),2);
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
        
        bubble_contours = tree_bubble_contours;
        for i=1:length(bubbles)
            if isempty(bubble_contours{i})
                indicators = (sum(dist(:,bubbles{i})<r(i),2)~=0) & (min(dist(:,bubbles{i})')'+r(i)/10<min(dist(:,setdiff(1:end,bubbles{i}))')');
                % it is possible that there are holes in the region labeled by indicators, we want to fill those holes so that the bubbles look nice.   
                tmp = reshape(double(indicators),sqrt(length(indicators)),sqrt(length(indicators)));
                tmp(1,tmp(1,:)==0)=-1;
                tmp(end,tmp(end,:)==0)=-1;
                tmp(tmp(:,1)==0,1)=-1;
                tmp(tmp(:,end)==0,end)=-1;
                old_tmp = tmp;
                while 1
                    tmp([-ones(1,size(tmp,2));tmp(1:end-1,:)]==-1 & tmp==0)=-1;
                    tmp([tmp(2:end,:);-ones(1,size(tmp,2))]==-1 & tmp==0)=-1;
                    tmp([-ones(size(tmp,1),1),tmp(:,1:end-1)]==-1 & tmp==0)=-1;
                    tmp([tmp(:,2:end),-ones(size(tmp,1),1)]==-1 & tmp==0)=-1;
                    if isequal(old_tmp,tmp)
                        break;
                    end
                    old_tmp=tmp;
                    % figure(1);subplot(2,2,4);imagesc(flipud(tmp))
                end
                indicators(tmp(:)==0)=1;
                
                outer_points = get_outer_points(grid_position, indicators, delta_x, delta_y);
                bubble_contours{i} = outer_points;
            end
        end
        tree_bubble_contours = bubble_contours;
        
        if length(tree_annotations)~=0
            for k=1:length(bubble_contours)
                if ~isempty(bubble_contours{k})
                    outer_points = bubble_contours{k};
                    line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
                end
            end
        end
        
    case 2
        'this option not available'
    otherwise
        'do nothing';
end


% % show node numbers
show_node_numbers = 0;
if show_node_numbers==1
    for k=1:size(coeff,2), text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7); end
end
axis off;


% % about colore
if show_color==0
    node_color = zeros(length(node_size),3) + 0.3;
elseif show_color==1
    node_color = zeros(length(node_size),3);
    color_code_data = node_data;
    % normalize color code data 
    % % this is unique to here 
    prc95 = prctile(data,98);
    prc05 = prctile(data,2); 
    color_code_data = (color_code_data - prc05)/(prc95-prc05);
    color_code_data(color_code_data<0.01)=0.01;
    color_code_data(color_code_data>0.99)=0.99;
    % define color map
    if ~exist('color_scheme')
        cmap_tmp = get_JET_color_map;
    elseif isequal(color_scheme,'jet')
        cmap_tmp = get_JET_color_map;
    elseif isequal(color_scheme,'gray')
        cmap_tmp = get_gray_color_map;
    elseif isequal(color_scheme,'half_jet')
        cmap_tmp = get_half_JET_color_map;
    else
        cmap_tmp = get_JET_color_map;
    end
    for k=1:size(coeff,2), 
        node_color(k,:) = interp1(((1:size(cmap_tmp,1))'-1)/(size(cmap_tmp,1)-1),cmap_tmp,color_code_data(k));  
        if sum(isnan(node_color(k,:)))~=0
            node_color(k,:) = [1,1,1];
        end
    end
else
    error('error in show_color value, it can be either 0 or 1');
end
    

% % draw nodes
for k=1:length(node_size)
    if exist('is_selected') && is_selected(k)==1
        % do nothing;
    else
        handle_tmp = plot(coeff(1,k),coeff(2,k),'o','markersize',node_size(k), 'markerfacecolor',node_color(k,:),'color',node_color(k,:));
    end
end
for k=1:length(node_size)
    if exist('is_selected') && is_selected(k)==1
        handle_tmp = plot(coeff(1,k),coeff(2,k),'s','markersize',node_size(k), 'markerfacecolor',node_color(k,:),'color',node_color(k,:));
    else
        % do nothing;
    end
end



% % control axis limits
axis_lims = reshape([-max(abs(coeff)'); +max(abs(coeff)')], 1, 4)*1.1;
for i=1:4
    if abs(axis_lims(i))<55
        if axis_lims(i)>=0
            axis_lims(i)=55;
        else
            axis_lims(i)=-55;
        end
    end
end
axis(axis_lims);



function SPADE_cmap = get_JET_color_map

SPADE_cmap = [   0         0    0.5625
                 0         0    0.6250
                 0         0    0.6875
                 0         0    0.7500
                 0         0    0.8125
                 0         0    0.8750
                 0         0    0.9375
                 0         0    1.0000
                 0    0.0625    1.0000
                 0    0.1250    1.0000
                 0    0.1875    1.0000
                 0    0.2500    1.0000
                 0    0.3125    1.0000
                 0    0.3750    1.0000
                 0    0.4375    1.0000
                 0    0.5000    1.0000
                 0    0.5625    1.0000
                 0    0.6250    1.0000
                 0    0.6875    1.0000
                 0    0.7500    1.0000
                 0    0.8125    1.0000
                 0    0.8750    1.0000
                 0    0.9375    1.0000
                 0    1.0000    1.0000
            0.0625    1.0000    0.9375
            0.1250    1.0000    0.8750
            0.1875    1.0000    0.8125
            0.2500    1.0000    0.7500
            0.3125    1.0000    0.6875
            0.3750    1.0000    0.6250
            0.4375    1.0000    0.5625
            0.5000    1.0000    0.5000
            0.5625    1.0000    0.4375
            0.6250    1.0000    0.3750
            0.6875    1.0000    0.3125
            0.7500    1.0000    0.2500
            0.8125    1.0000    0.1875
            0.8750    1.0000    0.1250
            0.9375    1.0000    0.0625
            1.0000    1.0000         0
            1.0000    0.9375         0
            1.0000    0.8750         0
            1.0000    0.8125         0
            1.0000    0.7500         0
            1.0000    0.6875         0
            1.0000    0.6250         0
            1.0000    0.5625         0
            1.0000    0.5000         0
            1.0000    0.4375         0
            1.0000    0.3750         0
            1.0000    0.3125         0
            1.0000    0.2500         0
            1.0000    0.1875         0
            1.0000    0.1250         0
            1.0000    0.0625         0
            1.0000         0         0
            0.9375         0         0
            0.8750         0         0
            0.8125         0         0
            0.7500         0         0
            0.6875         0         0
            0.6250         0         0
            0.5625         0         0
            0.5000         0         0];
        
    
function SPADE_cmap = get_half_JET_color_map

SPADE_cmap = [  0.5625    1.0000    0.4375
                0.6250    1.0000    0.3750
                0.6875    1.0000    0.3125
                0.7500    1.0000    0.2500
                0.8125    1.0000    0.1875
                0.8750    1.0000    0.1250
                0.9375    1.0000    0.0625
                1.0000    1.0000         0
                1.0000    0.9375         0
                1.0000    0.8750         0
                1.0000    0.8125         0
                1.0000    0.7500         0
                1.0000    0.6875         0
                1.0000    0.6250         0
                1.0000    0.5625         0
                1.0000    0.5000         0
                1.0000    0.4375         0
                1.0000    0.3750         0
                1.0000    0.3125         0
                1.0000    0.2500         0
                1.0000    0.1875         0
                1.0000    0.1250         0
                1.0000    0.0625         0
                1.0000         0         0
                0.9375         0         0
                0.8750         0         0
                0.8125         0         0
                0.7500         0         0
                0.6875         0         0
                0.6250         0         0
                0.5625         0         0
                0.5000         0         0];
        
    
function SPADE_cmap = get_gray_color_map
baseMap = [220   220    220;
           [165   165    165]*1.4;
           [110   110    110]*1.4;
           [55    55     55]*1.4;
           1     1      1]/255*(10/11);
idx1 = linspace(0,1,size(baseMap,1));
idx2 = linspace(0,1,100);
SPADE_cmap = interp1(idx1,baseMap,idx2);


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
%     if dummy > sqrt(delta_x^2+delta_y^2)*2;
%         flag_outer_points(I)=1; continue;
%     end
    ordered_points = [ordered_points; outer_points(I,:)]; flag_outer_points(I)=1;
end
outer_points = ordered_points;

