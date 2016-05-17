function [idx,c,c_ind]=deterministic_kmeans_for_SPADE(data,k,local_density)

fprintf('Deterministic k-means of the pooled downsampled cells ...');
tic;
[initial_idx, initial_Kmeans_centers] = initialKmeans(data,k);
[idx, c, c_ind]=kmeans_for_SPADE(data,local_density,k,initial_idx,initial_Kmeans_centers);
toc

end % function



function [idx,C, C_ind] = kmeans_for_SPADE(data,local_density,k,initial_idx,initial_Kmeans_centers)

maxit = 1000;

% ind_has_nan = find(sum(isnan(data),2)~=0);
% X = data;
% X(ind_has_nan,:)=[];

X = data;


% n points in p dimensional space
[n, p] = size(X);
% distance metric
distance = 'sqeuclidean';
C = initial_Kmeans_centers;
idx = initial_idx;
% D = distfun(X, C, distance);
% [d, idx] = min(D, [], 2);
m = accumarray(idx,1,[k,1]);
C_ind = zeros(max(idx),1);

moved = 1:n;            % which point moved to new cluster
changed = 1:k;          % which cluster is changed
previdx = zeros(n,1);   % previous idx
prevtotsumD = Inf;      % previous cost function total sum dist
iter = 0;
converged = false;
while true
    iter = iter + 1;

    % Calculate the new cluster centroids and counts, and update the
    % distance from every point to those new cluster centroids
    [C(changed,:), m(changed), C_ind(changed)] = gcentroids(X, idx, changed, distance, local_density);
    D(:,changed) = distfun(X, C(changed,:), distance);

    % Deal with clusters that have just lost all their members
    empties = changed(m(changed) == 0);
    if ~isempty(empties)
        error(message('stats:kmeans:EmptyCluster'));
        switch emptyact
            case 'drop'
                % Remove the empty cluster from any further processing
                D(:,empties) = NaN;
                changed = changed(m(changed) > 0);
            case 'singleton'
                for i = empties
                    d = D((idx-1)*n + (1:n)'); % use newly updated distances

                    % Find the point furthest away from its current cluster.
                    % Take that point out of its cluster and use it to create
                    % a new singleton cluster to replace the empty one.
                    [~, lonely] = max(d);
                    from = idx(lonely); % taking from this cluster
                    if m(from) < 2
                        % In the very unusual event that the cluster had only
                        % one member, pick any other non-singleton point.
                        from = find(m>1,1,'first');
                        lonely = find(idx==from,1,'first');
                    end
                    C(i,:) = X(lonely,:);
                    m(i) = 1;
                    idx(lonely) = i;
                    D(:,i) = distfun(X, C(i,:), distance, iter, rep, reps);

                    % Update clusters from which points are taken
                    [C(from,:), m(from)] = gcentroids(X, idx, from, distance, local_density);
                    D(:,from) = distfun(X, C(from,:), distance, iter, rep, reps);
                    changed = unique([changed from]);
                end
        end
    end

    % Compute the total sum of distances for the current configuration.
    totsumD = sum(D((idx-1)*n + (1:n)'));
    % Test for a cycle: if objective is not decreased, back out
    % the last step and move on to the single update phase
    if prevtotsumD <= totsumD
        idx = previdx;
        [C(changed,:), m(changed), C_ind(changed)] = gcentroids(X, idx, changed, distance, local_density);
        iter = iter - 1;
        break;
    end
    if iter >= maxit
        break;
    end

    % Determine closest cluster for each point and reassign points to clusters
    previdx = idx;
    prevtotsumD = totsumD;
    [d, nidx] = min(D, [], 2);

    % Determine which points moved
    moved = find(nidx ~= previdx);
    if ~isempty(moved)
        % Resolve ties in favor of not moving
        moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
    end
    if isempty(moved)
        converged = true;
        break;
    end
    idx(moved) = nidx(moved);

    % Find clusters that gained or lost members
    changed = unique([idx(moved); previdx(moved)])';

end % phase one



end % function


function [centroids, counts, centroids_ind] = gcentroids(X, index, clusts, dist, local_density)
%GCENTROIDS Centroids and counts stratified by group.
local_density = local_density(:)';
p = size(X,2);
num = length(clusts);
centroids = NaN(num,p);
counts = zeros(num,1);
centroids_ind = zeros(num,1);

for i = 1:num
    members = (index == clusts(i));
    if any(members)
        counts(i) = sum(members);
        switch dist
            case 'sqeuclidean'
                %centroids(i,:) = sum(X(members,:),1) / counts(i);
                tmp = sum(X(members,:).*repmat(local_density(members)',1,size(X,2)),1) / sum(local_density(members));
                [~,ind] = min(sum(abs(X(members,:)-repmat(tmp,counts(i),1)),2));
                centroids_ind(i) = max(find(members==1,ind));
                centroids(i,:) = X(centroids_ind(i),:);
            case 'cityblock'
                % Separate out sorted coords for points in i'th cluster,
                % and use to compute a fast median, component-wise
                members_local_density = local_density(members)';
                [Xsorted,I] = sort(X(members,:),1);
                [~,median_ind] = min(abs(cumsum(members_local_density(I))-sum(members_local_density)/2),[],1);
                centroids(i,:) = diag(Xsorted(median_ind,:));
        end
    end
end
end % function





function D = distfun(X, C, dist, iter,rep, reps)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);

switch dist
    case 'sqeuclidean'
        for i = 1:nclusts
            D(:,i) = (X(:,1) - C(i,1)).^2;
            for j = 2:p
                D(:,i) = D(:,i) + (X(:,j) - C(i,j)).^2;
            end
            % D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
        end
    case 'cityblock'
        for i = 1:nclusts
            D(:,i) = abs(X(:,1) - C(i,1));
            for j = 2:p
                D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
            end
            % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
        end
    case {'cosine','correlation'}
        % The points are normalized, centroids are not, so normalize them
        normC = sqrt(sum(C.^2, 2));
        if any(normC < eps(class(normC))) % small relative to unit-length data points
            if reps==1
                error(message('stats:kmeans:ZeroCentroid', iter));
            else
                error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
            end
            
        end
        
        for i = 1:nclusts
            D(:,i) = max(1 - X * (C(i,:)./normC(i))', 0);
        end
    case 'hamming'
        for i = 1:nclusts
            D(:,i) = abs(X(:,1) - C(i,1));
            for j = 2:p
                D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
            end
            D(:,i) = D(:,i) / p;
            % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
        end
end
end % function





function [idx, centers]=initialKmeans(data,k)
% method based on the following paper
% Su, Ting, and Jennifer Dy. "A deterministic method for initializing k-means clustering." Tools with Artificial Intelligence, 2004. ICTAI 2004. 16th IEEE International Conference on. IEEE, 2004.   

idx = ones(size(data,1),1);         % initialize the total number of group to be 1, every point in the same group
ssc = zeros(k,1);                   % allocate space for cluster SSC
centers = zeros(k,size(data,2));    % allocate space for cluster centers
current_num_groups = 1;             % initialize the total number of group to be 1
centers(current_num_groups,:) = mean(data(idx==current_num_groups,:),1);    
ssc(current_num_groups)=sum(sum(bsxfun(@minus,data(idx==current_num_groups,:),centers(current_num_groups,:)).^2));

while current_num_groups<k  % if the number of current groups is smaller than k, continue the loop and partition more clusters
    [~,idx_max_ssc_group]=max(ssc);                                         % find the index of the largest ssc group/cluster
    idx_points_in_max_ssc_group = find(idx==idx_max_ssc_group);             % find which points belong to the largest group/cluster
    
%     % divide points in the max ssc group by first principle component
%     [pc_eig_vec,~]=eigs(cov(data(idx_points_in_max_ssc_group,:)),1);        % find the principle component fo the largest scc group
%     if sum(pc_eig_vec)<0
%         pc_eig_vec=-pc_eig_vec;
%     end
%     threshold = centers(idx_max_ssc_group,:)*pc_eig_vec;                    % project the center to the principle component
%     subidx1=(data(idx_points_in_max_ssc_group,:)*pc_eig_vec<=threshold);    % project every point in the cluster to the principle component, and find out which points are to the left or right of the center after projection
    
    % divide points according to kmeans of the points along the first principle component 
    [pc_eig_vec,~]=eigs(cov(data(idx_points_in_max_ssc_group,:)),1);        % find the principle component fo the largest scc group
    if sum(pc_eig_vec)<0
        pc_eig_vec=-pc_eig_vec;
    end
    threshold = centers(idx_max_ssc_group,:)*pc_eig_vec;                    % project the center to the principle component
    subidx1_init = 1+(data(idx_points_in_max_ssc_group,:)*pc_eig_vec>threshold);    % project every point in the cluster to the principle component, and find out which points are to the left or right of the center after projection
    projected_data = data(idx_points_in_max_ssc_group,:)*pc_eig_vec;
    [subidx1,~,~] = kmeans_for_SPADE(projected_data,ones(1,length(projected_data)),2,subidx1_init,[min(projected_data);max(projected_data)]);
    subidx1 = (subidx1==1);

%     % divide points according to kmeans of the points along the first principle component 
%     [pc_eig_vec,~]=eigs(cov(data(idx_points_in_max_ssc_group,:)),1);        % find the principle component fo the largest scc group
%     if sum(pc_eig_vec)<0
%         pc_eig_vec=-pc_eig_vec;
%     end
%     projected_data = data(idx_points_in_max_ssc_group,:)*pc_eig_vec;
%     obj = gmdistribution.fit(projected_data,2);   
%     subidx1 = (cluster(obj,projected_data)==1);


    % subplot(1,2,1); data2 = data'; for i=1:max(idx), plot(data2(1,idx==i),data2(2,idx==i),'.','color',[rand,rand,rand]); tmp =  [data2(1,idx==i);data2(2,idx==i)]; h = convhull(data2(1,idx==i),data2(2,idx==i)); line(tmp(1,h),tmp(2,h)); hold on; end
    idx(idx_points_in_max_ssc_group(subidx1)) = current_num_groups + 1;     % creat a new cluster
    % subplot(1,2,2); data2 = data'; for i=1:max(idx), plot(data2(1,idx==i),data2(2,idx==i),'.','color',[rand,rand,rand]); tmp =  [data2(1,idx==i);data2(2,idx==i)]; h = convhull(data2(1,idx==i),data2(2,idx==i)); line(tmp(1,h),tmp(2,h)); hold on; end
    centers(idx_max_ssc_group,:) = mean(data(idx==idx_max_ssc_group,:),1);  % computer new center and ssc for the one old and one new cluster after the partitioning
    ssc(idx_max_ssc_group)=sum(sum(bsxfun(@minus,data(idx==idx_max_ssc_group,:),centers(idx_max_ssc_group,:)).^2));
    centers(current_num_groups + 1,:) = mean(data(idx==current_num_groups + 1,:),1);    
    ssc(current_num_groups + 1)=sum(sum(bsxfun(@minus,data(idx==current_num_groups + 1,:),centers(current_num_groups + 1,:)).^2));
    current_num_groups = current_num_groups + 1;                            % update the total number of clusters
end
end % function


