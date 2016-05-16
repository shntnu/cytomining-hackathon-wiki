function [sample_group_progression_tree, idx_samples, sample_modules] = SPADE_cluster_cells(data, target_num_nodes)

% this is the same as group_flow_samples_to_tree.m in the older versions

idx_samples = 1:size(data,2);
counter = 1;
while max(idx_samples)>=1.5*target_num_nodes
    fold = max(2,round(max(idx_samples)/5000)); 
    fprintf('merging %d groups into about %d groups, round %d\n', max(idx_samples), ceil(max(idx_samples)/fold), counter);
    tmp_ind = find(idx_samples~=0);
    if isempty(find(idx_samples==0))
        idx_samples = merge_sample_groups_old(data, fold,idx_samples);
    else
        idx_samples(tmp_ind) = merge_sample_groups_old(data(:,tmp_ind), fold,idx_samples(tmp_ind));
    end
    if counter==5 % after 5 rounds of iterations, if a node still contains only one cell, this node is thrown away
        for i=1:max(idx_samples)
            if sum(idx_samples==i)==1
                idx_samples(idx_samples==i)=0;
            end
        end
        idx_samples = standardize_idx(idx_samples);
        display(['throw away ',num2str(sum(idx_samples==0)),' singltons'])
    end
    counter = counter + 1;
end 

sample_modules = get_module_mean(data',idx_samples)';
sample_group_progression_tree = mst(sample_modules');

return



function idx_new = merge_sample_groups_old(data, fold, idx)
% the samples grouped into give N groups by idx, 
% reduce the number of groups by a factor of fold (min of fold is 2)
% look at the groups one by one, 
% for each group, use single linkage to find the closet other group
% merge them together, and this other group is out of the game

isactive_module = ones(1,max(idx));
module_size = zeros(1,max(idx)); for i=1:length(module_size), module_size(i) = sum(idx==i); end
isactive_sample = ones(1,size(data,2));
fprintf('merging in progress ... groups left: %7d',max(idx));
iter = 1; total_num_modules = max(idx);
while sum(isactive_module)>1
%     module_ind = randsample(find(isactive_module==1),1);
% tic  
    tmp = module_size; tmp(module_size==0)=Inf;  tmp(isactive_module==0) = Inf;
    candidate_module_ind = find(tmp==min(tmp(isactive_module==1)) & isactive_module==1);
    if length(candidate_module_ind)==1
        module_ind = candidate_module_ind;
    else
        module_ind = candidate_module_ind(1); %randsample(candidate_module_ind,1);
    end

    sample_in_module = find(idx==module_ind);
    module_center = median(data(:,sample_in_module),2);
    dist = sum(abs(repmat(module_center,1,size(data,2)) - data),1);
    dist(idx==module_ind) = Inf;    % dist to everyone in my own group
    [Y,I] = sort(dist,'ascend'); 
    module_to_be_deleted = idx(I(1:fold-1));
    first_inactive_module_ind = find(isactive_module(module_to_be_deleted)==0,1);
    if ~isempty(first_inactive_module_ind) && first_inactive_module_ind==1, module_to_be_deleted=[]; end
    if ~isempty(first_inactive_module_ind) && first_inactive_module_ind~=1, module_to_be_deleted=module_to_be_deleted(1:first_inactive_module_ind-1); module_to_be_deleted = unique(module_to_be_deleted); end
% toc
% fprintf('%7d %7d ',total_num_modules,sum(isactive_module==1));drawnow;
    if isempty(module_to_be_deleted)
        isactive_module(module_ind)=0;
%         if sum(isactive_module==1)<2000
%             1;
%         end
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d %7d ',total_num_modules,sum(isactive_module==1));drawnow;
        continue;
    end
    isactive_module(module_to_be_deleted)=0;
    isactive_module(module_ind)=0;
    for i=1:length(module_to_be_deleted)
        idx(idx==module_to_be_deleted(i)) = module_ind;  % merge
    end
    module_size(module_to_be_deleted)=0; module_size(module_ind) = sum(idx==module_ind);
    isactive_sample(idx==module_ind)=0;
    total_num_modules = total_num_modules - length(module_to_be_deleted);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d %7d ',total_num_modules,sum(isactive_module==1));
    iter = iter + 1; % iter is only serving for the drawnow in the following lines
    if round(iter/10)==iter/10
        drawnow;
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d %7d',total_num_modules,sum(isactive_module==1));drawnow;
fprintf('\n')
idx_new = standardize_idx(idx);
return    


function idx_new = merge_sample_groups_new(data, fold, idx)
% the samples grouped into give N groups by idx, 
% reduce the number of groups by a factor of fold (min of fold is 2)
% look at the groups one by one, 
% for each group, use single linkage to find the closet other group
% merge them together, and this other group is out of the game

isactive_module = ones(1,max(idx));
isactive_sample = ones(1,size(data,2));
fprintf('merging in progress ... groups left: %7d',max(idx));
iter = 1; total_num_modules = max(idx);
while sum(isactive_module)>1
    module_ind = randsample(find(isactive_module==1),1);
    sample_in_module = find(idx==module_ind);
    module_center = median(data(:,sample_in_module),2);
    dist = zeros(1,length(isactive_module))+Inf;
    for i=1:length(dist)
        if i==module_ind, continue; end
        if isactive_module(i)==0, continue; end
        dist(i) = sum(abs(module_center - median(data(:,idx==i),2)));
    end
    [Y,I] = sort(dist,'ascend'); 
    module_to_be_deleted = I(1:fold-1);
    isactive_module(module_to_be_deleted)=0;
    isactive_module(module_ind)=0;
    for i=1:length(module_to_be_deleted)
        idx(idx==module_to_be_deleted(i)) = module_ind;  % merge
    end
    isactive_sample(idx==module_ind)=0;
    total_num_modules = total_num_modules - length(module_to_be_deleted);
    fprintf('\b\b\b\b\b\b\b%7d',total_num_modules);
    iter = iter + 1; % iter is only serving for the drawnow in the following lines
    if round(iter/50)==iter/50
        drawnow;
    end
end
fprintf('\n')
idx_new = standardize_idx(idx);
return    