function [center_dist_scaled, center_dist, scaling_factor]= SPADE_center_pairwise_dist(idx,centers,centers_ind, data_points, local_density) 
% each column is one point

center_dist = squareform(pdist(centers'));
center_dist_scaled = center_dist;
scaling_factor = zeros(size(center_dist));

tic
fprintf('compute scaled cluster distance by density ... %3d%%',0);
for i=1:size(center_dist,1)-1
    for j=i+1:size(center_dist,1)
        v1 = data_points(:,centers_ind(i));
        v2 = data_points(:,centers_ind(j));
        projection_to_segment = (v2'-v1')*((repmat(v2,1,size(data_points,2))-data_points))/((v2'-v1')*(v2-v1));
        projection_squared_error = sum( (v1*projection_to_segment + v2*(1-projection_to_segment) - data_points).^2, 1);
        [f,xi] = ksdensity(projection_to_segment( projection_squared_error < (norm(v2-v1)/2)^2 ),[0:0.01:1]);
        center_dist_scaled(i,j) = center_dist(i,j) * (1+ log(min(f([1,end])) / (min(f)+eps)));
        center_dist_scaled(j,i) = center_dist_scaled(i,j);
        scaling_factor(i,j) = (1+ log(min(f([1,end])) / (min(f)+eps)));
        scaling_factor(j,i) = scaling_factor(i,j);
    end
    fprintf('\b\b\b\b%3d%%',round(100*i/(size(center_dist,1)-1)));
end
toc


% tic
% fprintf('compute scaled cluster distance by density ... %3d%%',0);
% for i=1:size(center_dist,1)-1
%     for j=i+1:size(center_dist,1)
%         ind = (idx==i | idx==j);
%         v1 = data_points(:,centers_ind(i));
%         v2 = data_points(:,centers_ind(j));
%         data_tmp = data_points(:,ind);
%         projection_to_segment = (v2'-v1')*((repmat(v2,1,size(data_tmp,2))-data_tmp))/((v2'-v1')*(v2-v1));
%         [f,xi] = ksdensity(projection_to_segment,[0:0.01:1]);
%         center_dist_scaled(i,j) = center_dist(i,j) * (1+ log(min(f([1,end])) / (min(f)+eps)));
%         center_dist_scaled(j,i) = center_dist_scaled(i,j);
%     end
%     fprintf('\b\b\b\b%3d%%',round(100*i/(size(center_dist,1)-1)));
% end
% toc