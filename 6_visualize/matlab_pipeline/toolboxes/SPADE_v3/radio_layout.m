function positions = radio_layout(mst_tree, data)

% make sure that tree_adj is undirected
adj = triu(mst_tree,1); adj = adj + adj';


% find the backbone, which is the longest path
shortest_hop = tree_shortest_hop(adj);
shortest_hop = shortest_hop - diag(diag(shortest_hop));
[ind_j,ind_i]  = find_matrix_top_element(shortest_hop);
ind_i = ind_i(1); ind_j = ind_j(1);
back_bones = find(shortest_hop(ind_i,:)+shortest_hop(ind_j,:)==shortest_hop(ind_i,ind_j));
[dummy,I] = sort(shortest_hop(ind_i,back_bones));
back_bones = back_bones(I);


% create parent_list representation of the tree using middle of backbone as root
root = back_bones(round(length(back_bones)/2));
N = size(adj,1);
parent_list = -1 + zeros(N,1);
depth = -1 + zeros(N,1);
parent_list(root)=0;
depth(root)=0;
queue = root;
while ~isempty(queue)
    current_node = queue(1);
    reachable_unlabeled_nodes = find(adj(:,current_node)==1 & (parent_list == -1));
    parent_list(reachable_unlabeled_nodes) = current_node;
    depth(reachable_unlabeled_nodes) = depth(current_node)+1;
    queue = [queue;reachable_unlabeled_nodes(:)];
    queue(1) = [];
end


% % convert parent list to undirected graph adj matrix
% adj2 = zeros(N,N);
% edge_list = [];
% for i=1:length(parent_list)
%     if parent_list(i)==0
%         continue;
%     else
%         adj2(i, parent_list(i))=1;
%         adj2(parent_list(i), i)=1;
%         edge_list = [edge_list; [parent_list(i), i]];
%     end
% end

fprintf('  Compute initial radio expansion layout ... \n');

dist = pdist(data','cityblock'); clustTreeEuc = linkage(dist,'single');
% dist = pdist(data','euclidean'); clustTreeEuc = linkage(dist,'average');
h=figure(107); [H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2)); close(h);

% inspired by matlab treeplot, radioplot.  Although there is no gaurantee that there will be no edge crossing, it is unlikely.  
node_data = zeros(size(data));
root = find(parent_list==0);
queue = root;
angle_center = zeros(1,N); angle_center(root) = 0;
angle_range = zeros(1,N); angle_range(root)  = 2*pi;
node_data(:,root) = mean(data(:,:),2);
x = zeros(1,N);
y = zeros(1,N);
while ~isempty(queue)
    current_node = queue(1);
    node_data(:,current_node) = data(:,current_node);
    reachable_nodes = find(parent_list==current_node);
    reachable_nodes = perm_c(ismember(perm_c, reachable_nodes));
    num_leaf_nodes = zeros(length(reachable_nodes),1);
    num_depth = zeros(length(reachable_nodes),1);
    num_subtree_nodes = zeros(length(reachable_nodes),1);
    for i=1:length(reachable_nodes)
        [num_leaf_nodes(i), num_depth(i), num_subtree_nodes(i), subtree_nodes] = find_num_leaf_nodes(parent_list, reachable_nodes(i));
        node_data(:,reachable_nodes(i)) = mean(data(:,subtree_nodes),2);
    end
    new_angle_range =  angle_range(current_node)/sum(num_leaf_nodes)*num_leaf_nodes;
    tmp = cumsum(reshape([new_angle_range(:)'/2; new_angle_range(:)'/2],1,length(reachable_nodes)*2));
    new_angle_center = angle_center(current_node)-angle_range(current_node)/2 + tmp(1:2:end);
    for i=1:length(reachable_nodes)
        x(reachable_nodes(i)) = x(current_node) + cos(new_angle_center(i));
        y(reachable_nodes(i)) = y(current_node) + sin(new_angle_center(i));
    end
    queue = [queue;reachable_nodes(:)];
    angle_range(reachable_nodes) = new_angle_range(:);
    angle_center(reachable_nodes) = new_angle_center(:);
    queue(1) = [];
end
% subplot(1,2,1); display_graph(adj,[x;y])

positions = [x;y];



% rotation to make the tree more spreadout
iter = 0;
fprintf('  Refine the layout by rotation of subtrees (max 1000 iter) ... %5d',iter);
while 1

    all_rotation_forces = [];
    for i=1:N
        if i==root
            continue; 
        end
        current_node = i;
        rotation_force = compute_rotation_force(current_node, parent_list, positions);
        all_rotation_forces(i) = rotation_force;
    end
    
    for i=1:N
        [max_force_amplitude,current_node] = max(abs(all_rotation_forces));
        if max_force_amplitude==0
            break;
        end
        subtree_nodes = find_subtree_nodes(parent_list, current_node);  % find out all the nodes that are under the current node, including itself
        % rotate pi/180 and see whether things become better (force become smaller) 
        new_positions = positions;
        if all_rotation_forces(current_node)>0
            amout_rotation = pi/180*1;
        else
            amout_rotation = -pi/180*1;
        end
        new_positions(:,subtree_nodes) =  [cos(amout_rotation) sin(amout_rotation); -sin(amout_rotation) cos(amout_rotation)] * (new_positions(:,subtree_nodes) - repmat(positions(:,parent_list(current_node)),1,length(subtree_nodes)))  + repmat(positions(:,parent_list(current_node)),1,length(subtree_nodes));
        new_rotation_force = compute_rotation_force(current_node, parent_list, new_positions);
        angle_tuned = 0;
        if abs(new_rotation_force)<abs(all_rotation_forces(current_node))
            while abs(new_rotation_force)<abs(all_rotation_forces(current_node))
                positions = new_positions;
                all_rotation_forces(current_node) = new_rotation_force;
                new_positions(:,subtree_nodes) =  [cos(amout_rotation) sin(amout_rotation); -sin(amout_rotation) cos(amout_rotation)] * (new_positions(:,subtree_nodes) - repmat(positions(:,parent_list(current_node)),1,length(subtree_nodes)))  + repmat(positions(:,parent_list(current_node)),1,length(subtree_nodes));
                new_rotation_force = compute_rotation_force(current_node, parent_list, new_positions);
                angle_tuned = angle_tuned + 1;
            end
        end
        if angle_tuned >=2 
            break;
        else
            all_rotation_forces(current_node)=0;
        end
    end
    
%     subplot(1,2,2); display_graph(adj,new_positions)
%     text(new_positions(1,current_node), new_positions(2,current_node), num2str(current_node))
%     drawnow;
    
    if max_force_amplitude==0
        break;
    end
    
    iter = iter+1;
    fprintf('\b\b\b\b\b%5d',iter);
    if iter>=1000 
        break;
    end
end





function rotation_force = compute_rotation_force(current_node, parent_list, positions)
% tic;
N = length(parent_list);
subtree_nodes = find_subtree_nodes(parent_list, current_node);  % find out all the nodes that are under the current node, including itself
outside_nodes = setdiff(1:N, subtree_nodes);

rotation_force = 0;
for j = 1:length(subtree_nodes)  % for each node under the current subtree, figure out the force and rotation force
    % figure out the repelling force on this onde
    tmp_diff = repmat(positions(:,subtree_nodes(j)),1,length(outside_nodes))-positions(:,outside_nodes);
    tmp_norm = sqrt(sum(tmp_diff.^2));
    tmp_magnitude =  1./(tmp_norm.^2); tmp_magnitude(tmp_norm>10)=0;
    force = sum(tmp_diff.*repmat(tmp_magnitude./tmp_norm,2,1),2);
    % translate the repelling force in to rotation force
    radius_vector = positions(:,subtree_nodes(j)) - positions(:,parent_list(current_node));
    radius_length = norm(radius_vector);
    radius_direction = radius_vector/radius_length;
    orthoganal_force = force - (force'*radius_direction)*radius_direction;
    if ([cos(pi/2) sin(pi/2); -sin(pi/2) cos(pi/2)]*radius_direction)'*orthoganal_force > 0 % clockwise
        rotation_force = rotation_force + norm(orthoganal_force)*radius_length;
    else % counter-clockwise
        rotation_force = rotation_force - norm(orthoganal_force)*radius_length;
    end
end
% toc


function [num_leaf_nodes, num_depth, num_subtree_nodes, subtree_nodes] = find_num_leaf_nodes(parent_list, start_node)
queue = start_node;
depth = 1;
num_leaf_nodes = 0;
num_depth = 1;
num_subtree_nodes = 1;
subtree_nodes = start_node;
while ~isempty(queue)
    current_node = queue(1);
    reachable_nodes = find(parent_list==current_node);
    subtree_nodes = [subtree_nodes; reachable_nodes(:)];
    num_subtree_nodes = num_subtree_nodes + length(reachable_nodes(:));
    if isempty(reachable_nodes)
        num_leaf_nodes = num_leaf_nodes + 1;
    else
        queue = [queue;reachable_nodes(:)];
        depth = [depth;repmat(depth(1)+1, length(reachable_nodes),1)];
    end
    num_depth = max(union(num_depth, depth));
    depth(1)=[];
    queue(1)=[];
end


function subtree_nodes = find_subtree_nodes(parent_list, start_node)
queue = start_node;
subtree_nodes = start_node;
while ~isempty(queue)
    current_node = queue(1);
    reachable_nodes = find(parent_list==current_node);
    subtree_nodes = [subtree_nodes; reachable_nodes(:)];
    queue = [queue;reachable_nodes(:)];
    queue(1)=[];
end





function I = order_reachable_nodes_to_circle(node_data)

N = size(node_data,2);
dist = pdist(node_data','cityblock');
dist = squareform(dist);
for i=1:N-1  % N-1 edges, just one before forming circle
    
end

    

