function [ind_i,ind_j] = find_matrix_top_element(adj)

ind = find(adj==max(max(adj)));
dim = size(adj,1);
for i=1:length(ind)
    ind_j(i) = ceil(ind(i)/dim);
    ind_i(i) = ind(i) - (ind_j(i)-1)*dim;
end
return
