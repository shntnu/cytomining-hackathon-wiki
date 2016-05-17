function regions = get_all_possible_tree_regions(mst_tree)
% regions = get_all_possible_tree_regions(mst_tree)
% NOTE: output is regions by nodes, each row is an indicator vector, where
%       1 means this node belongs to this region

[pairs] = find_matrix_big_element(triu(abs(mst_tree),1), 1e-5);
regions = zeros(size(mst_tree,1), size(pairs,1)*2+size(pairs,1)*(size(pairs,1)-1)/2*3);
count = 1;
for i=1:size(pairs,1) % remove 1 edge, and look at pieces after cut
    adj = mst_tree;
    adj(pairs(i,1),pairs(i,2))=0;
    adj(pairs(i,2),pairs(i,1))=0;
    [regulation_matrix, kept_ind, components_size,components] = extract_connected_component(adj);
    regions(:,count:count+1) = components;
    count = count + 2;
end
for i=1:size(pairs,1)-1 % remove 2 edges, and look at pieces after cut
    for j=i+1:size(pairs,1)
        adj = mst_tree;
        adj(pairs(i,1),pairs(i,2))=0;
        adj(pairs(i,2),pairs(i,1))=0;
        adj(pairs(j,1),pairs(j,2))=0;
        adj(pairs(j,2),pairs(j,1))=0;
        [regulation_matrix, kept_ind, components_size,components] = extract_connected_component(adj);
        regions(:,count:count+2) = components;
        count = count + 3;
    end
end

regions = unique(regions','rows');
return


