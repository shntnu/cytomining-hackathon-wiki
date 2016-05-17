function draw_SPADE_tree(mst_tree, node_positions, node_size, node_data, data, show_color, show_annotation,is_selected, color_scheme)
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
        'show all bubbles';
    case 2 
        'show one selected bubble';
    otherwise
        3;
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
   