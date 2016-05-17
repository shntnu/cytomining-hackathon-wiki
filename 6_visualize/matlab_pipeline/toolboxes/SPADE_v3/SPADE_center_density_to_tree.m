function [A]= SPADE_center_density_to_tree(idx,centers,centers_ind, data_points, local_density) 
% each column is one point

center_dist = squareform(pdist(centers','cityblock'));  % center_dist = squareform(pdist(centers'));
A = mst_from_dist_matrix(center_dist); 
k=3;
for i=1:size(A,1)
    [~,ind] = sort(center_dist(i,:),'ascend');
    for j=2:k+1
        A(i,ind(j))=1;
        A(ind(j),i)=1;
    end
end

center_dist_scaled = center_dist.*A;

tic
fprintf('compute scaled cluster distance by density ... %3d%%',0);
for i=1:size(center_dist,1)-1
    for j=i+1:size(center_dist,1)
        if center_dist_scaled(i,j)~=0
            v1 = data_points(:,centers_ind(i));
            v2 = data_points(:,centers_ind(j));
            projection_to_segment = (v2'-v1')*((repmat(v2,1,size(data_points,2))-data_points))/((v2'-v1')*(v2-v1));
            projection_squared_error = sum( abs(v1*projection_to_segment + v2*(1-projection_to_segment) - data_points), 1);
            [f,xi] = ksdensity(projection_to_segment( projection_squared_error < sum(abs(v2-v1))/2 ),[0:0.01:1]);
            (1+ log(min(f([1,end])) / (min(f)+eps)))
            center_dist_scaled(i,j) = center_dist(i,j) * (1+ log(min(f([1,end])) / (min(f)+eps)));
            center_dist_scaled(j,i) = center_dist_scaled(i,j);
        end
    end
    fprintf('\b\b\b\b%3d%%',round(100*i/(size(center_dist,1)-1)));
end
toc

center_dist_scaled(center_dist_scaled==0) = max(center_dist_scaled(:))*100;
A = mst_from_dist_matrix(center_dist_scaled); 



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