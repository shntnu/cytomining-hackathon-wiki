function [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(individual_score, idx)
% individual_score and idx are both vectors of the same length

individual_score=individual_score(:);
idx = idx(:);

[dummy,reordering] = sort(idx);
idx = idx(reordering);
individual_score = individual_score(reordering);

group_seg_ends = [find(idx(1:end-1)~=idx(2:end)); length(idx)];
group_seg_starts = [1; find(idx(1:end-1)~=idx(2:end))+1];
counts = group_seg_ends-[0; group_seg_ends(1:end-1)];
group_avg = zeros(1,length(counts));
group_idx_values = zeros(1,length(counts));
for i=1:length(group_seg_ends)
    group_avg(i) = nanmean(individual_score(group_seg_starts(i):group_seg_ends(i)));
    group_idx_values(i) = unique(idx(group_seg_starts(i):group_seg_ends(i)));
end



% % the following code was created when I was looking at next-gen seq data,
% % where the individual scores are not NaNs
% 
% cumsum_scores = cumsum(individual_score);
% group_seg_ends = [find(idx(1:end-1)~=idx(2:end)); length(idx)];
% 
% counts = group_seg_ends-[0; group_seg_ends(1:end-1)];
% sum_scores = cumsum_scores(group_seg_ends) - [0; cumsum_scores(group_seg_ends(1:end-1))];
% 
% group_avg = sum_scores./counts;


