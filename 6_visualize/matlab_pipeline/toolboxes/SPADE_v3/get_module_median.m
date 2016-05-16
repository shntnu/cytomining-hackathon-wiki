function [module_median] = get_module_median(data, idx)
% [module_median] = get_module_median(data, idx)
% [module_median] = get_module_median(data, components)
% idx or components are outputs of cut_corr_matrix.m

if size(idx,2)~=1 && size(idx,1)~=1
    % idx is actually components
    module_median = zeros(size(idx,2),size(data,2));
    for i=1:size(idx,2)
        ind = find(idx(:,i)==1);
        module_median(i,:) = median(data(ind,:),1);
    end
        
else
    possible_clusters = setdiff(unique(idx),0);
    module_median = zeros(length(possible_clusters),size(data,2));
    for i=1:length(possible_clusters)
        ind = find(idx==possible_clusters(i));
        module_median(i,:) = median(data(ind,:),1);
    end
    % module_data = per_gene_normalization(module_data);
end
