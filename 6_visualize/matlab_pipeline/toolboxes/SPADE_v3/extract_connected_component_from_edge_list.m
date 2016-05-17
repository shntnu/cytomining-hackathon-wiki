function [edge_list, kept_ind, components_size,components] = extract_connected_component_from_edge_list(edge_list)
% [regulation_matrix, kept_ind, components_size,components] = extract_connected_component_from_edge_list(edge_list)

N = max(edge_list(:));
edge_list = unique([edge_list;repmat((1:N)',1,2)], 'rows');  % add self edges which is useful  

components=[];
is_assigned = zeros(N,1);
while sum(is_assigned==0)~=0
%     seed_nodes = find(is_assigned==0,1);
%     reachable_nodes = unique(edge_list(ismember(edge_list(:,1),seed_nodes) | ismember(edge_list(:,2),seed_nodes),:));
%     while ~isequal(seed_nodes,reachable_nodes)
%         seed_nodes = reachable_nodes;
%         reachable_nodes = unique(edge_list(ismember(edge_list(:,1),seed_nodes) | ismember(edge_list(:,2),seed_nodes),:));
%     end
%     reachable_nodes_old = seed_nodes;
    seed_nodes = find(is_assigned==0,1);
    reachable_new_nodes = setdiff(unique(edge_list(ismember(edge_list(:,1),seed_nodes) | ismember(edge_list(:,2),seed_nodes),:)), seed_nodes);
    while ~isempty(reachable_new_nodes)
        seed_nodes = union(seed_nodes, reachable_new_nodes);
        reachable_new_nodes = setdiff(unique(edge_list(ismember(edge_list(:,1),reachable_new_nodes) | ismember(edge_list(:,2),reachable_new_nodes),:)), seed_nodes);
    end
    e = zeros(N,1);
    e(seed_nodes)=1;
    components = [components,e];
    is_assigned(seed_nodes)=1;
end

components_size = sum(components);
component_ind = find(components_size == max(components_size),1);

kept_ind = find(components(:,component_ind)==1);
edge_list = edge_list(ismember(edge_list(:,1),kept_ind) & ismember(edge_list(:,2),kept_ind) , :);

return