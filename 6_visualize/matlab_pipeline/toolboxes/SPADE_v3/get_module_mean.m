function [module_mean] = get_module_mean(data, idx)
% [module_mean] = get_module_mean(data, idx)
% [module_mean] = get_module_mean(data, components)
% idx or components are outputs of cut_corr_matrix.m

if size(idx,2)~=1 && size(idx,1)~=1
    % idx is actually components
    module_mean = zeros(size(idx,2),size(data,2));
    for i=1:size(idx,2)
        ind = find(idx(:,i)==1);
        module_mean(i,:) = mean(data(ind,:),1);
    end
        
else
    possible_clusters = setdiff(unique(idx),0);
    module_mean = zeros(length(possible_clusters),size(data,2));

    for i=1:length(possible_clusters)
        ind = find(idx==possible_clusters(i));
        module_mean(i,:) = mean(data(ind,:),1);
    end
    % module_data = per_gene_normalization(module_data);
end
