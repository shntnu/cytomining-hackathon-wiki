function [is_connected, shortest_hop] = check_graph_connectivity(adj_matrix)
% [is_connected, shortest_hop] = check_graph_connectivity(adj_matrix)
% for an given adj_matrix of a undirected graph, first calculate the 
% shortest_hop between each pair of nodes. If any entry of the matrix 
% shortest_hop is zero, meaning that there is no connection between this 
% pair of nodes, we know that the graph is not connected, then this program
% output the number of subgraphs, and the number of nodes in each subgraph

adj_matrix = (abs(adj_matrix) + abs(adj_matrix') + eye(size(adj_matrix))) >0;
nNodes=size(adj_matrix,1);

shortest_hop = zeros(nNodes);
for i=1:nNodes
    e=zeros(nNodes,1); e(i)=1; hop=0;
    while sum(shortest_hop(i,:)~=0)<nNodes
        e = double(adj_matrix*e >0); hop=hop+1;
        ind = find(e'>0 & shortest_hop(i,:)==0);
        if isempty(ind), break; end
        shortest_hop(i,ind) = hop;
    end
end


% if not connected, calculate the size of each component
if sum(shortest_hop(1,:)~=0)~=size(shortest_hop,2)
    is_connected=0;
    component_size=[];
    count=0;
    tmp = shortest_hop; rearrange_shortest_hop=[];
    while sum(tmp(1,:)==0)>0
        ind = find(tmp(1,:)~=0); count=count+1;component_size=[component_size, length(ind)];
        rearrange_shortest_hop = [rearrange_shortest_hop, zeros(size(rearrange_shortest_hop,1),length(ind));zeros(length(ind),size(rearrange_shortest_hop,1)),tmp(ind,ind)];
        tmp=tmp(setdiff(1:size(tmp,1),ind),setdiff(1:size(tmp,1),ind));
    end
    count=count+1;
    rearrange_shortest_hop = [rearrange_shortest_hop, zeros(size(rearrange_shortest_hop,1),size(tmp,2));zeros(size(tmp,1),size(rearrange_shortest_hop,1)),tmp];
    component_size=[component_size, size(tmp,1)];
    count
    fliplr(sort(component_size))
else
    is_connected=1;
end
    return
