function target_density = downsample_to_certain_num_cells(data, local_density, desired_num)
% keep_prob = x./local_density
% need to find the value of "x", such that if we downsample according to
% "keep_prob", we end up with about "desired_num" cells
% therefore, need to solve the following
%      sum(min(x/local_density(i),1)) = desired_num
% which is equivalent to
%      x = (desired_num-i) / sum(1/local_density(i+1:end)) && local_density(i)<=x<=local_density(i+1) 

if desired_num>=length(local_density)
    target_density = max(local_density)+1;
    return
end
ld = [sort(local_density,'ascend')];
if desired_num/sum(1./local_density) <= ld(1)
    target_density = desired_num/sum(1./local_density);
    return
end
for i=1:length(ld)-1
    x = (desired_num-i) / sum(1./ld(i+1:end));
    if ld(i)<=x && x<=ld(i+1) 
        break;
    end
end
target_density = x;
return