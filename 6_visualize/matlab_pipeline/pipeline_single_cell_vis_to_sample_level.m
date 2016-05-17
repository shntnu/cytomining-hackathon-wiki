%% addpath to the matlab functions needed for this pipeline
addpath(genpath('toolboxes'))


%% read data source folder and get all file names
data_folder = 'single_cell_data';
filenames = getfilenames(data_folder);
filenames = filenames(isInListEnd(filenames,'.txt'));




%% read data, downsample, pooling
pooled_data = [];
pooled_downsampled_data = [];
MOA = []; % perturbation Mechanism Of Action
for i=1:length(filenames)
    MOA = [MOA; {filenames{i}(1:find(filenames{i}=='-',1)-1)}];            % in this example, MOA of each sample is encoded in the filenames
    tmp = read_xls_v2(fullfile(data_folder,filenames{i}), char(9),0);      % this read file function is obtained from the matlab implementation of SPADE
    feature_names = tmp(1,:)';                                             % in txt file for each sample, first row contains feature names, 
    tmp_data = cell2num(tmp(2:end,:));                                     % each subsequent row is a cell, each column is a feature
    pooled_data = [pooled_data; [tmp_data,ones(size(tmp_data,1),1)*i]];    % pool all data files together
    if size(tmp_data,1)>1000                                               % optional if data size too large, we many want to downsample the data     
        ind = sort(randsample(1:size(tmp_data,1),1000));                   % here is an example for simple uniform downsampling. More sophisticated methods such as density dependent downsampling or faithfull downsamples can be found in SPADE and SamSpectral
        tmp_data = tmp_data(ind,:);                                        % NOTE: if downsampleing is performed, we need to upsample in all subsequent visualizations. How to upsample will be different for each visualization algorithm. Examples are not provided here.
        pooled_downsampled_data = [pooled_downsampled_data; [tmp_data,ones(size(tmp_data,1),1)*i]];
    end
end


%% PCA-2D
figure(1)
data = pooled_data(:,1:end-1);
labels = pooled_data(:,end);
mappedX = multiclass_pcaplot(data',[],[1 2])'; % run PCA implemented in the SPADE package
subplot(1,2,1); plot(mappedX(:,1),mappedX(:,2),'g.')
subplot(1,2,2); FlowJo_contour2D(mappedX(:,1),mappedX(:,2),1)
figure(2)
counter = 1;
for i=[8 11 14 62 63 65]  % visualize 6 selected samples from the total of 103 samples. These 6 belong two MOAs. The first three are Aurora Kinase inhibitors, and the last three are labeled as general Kinase inhibitors.
    subplot(2,3,counter); counter = counter+1;
    plot(mappedX(:,1),mappedX(:,2),'g.',mappedX(labels==i,1),mappedX(labels==i,2),'b.');
    
    tmp1 = prctile(mappedX(:,1),1)  - (prctile(mappedX(:,1),99) - prctile(mappedX(:,1),1))/3;
    tmp2 = prctile(mappedX(:,1),99) + (prctile(mappedX(:,1),99) - prctile(mappedX(:,1),1))/3;
    xlim([tmp1, tmp2])
    tmp1 = prctile(mappedX(:,2),1)  - (prctile(mappedX(:,2),99) - prctile(mappedX(:,2),1))/3;
    tmp2 = prctile(mappedX(:,2),99) + (prctile(mappedX(:,2),99) - prctile(mappedX(:,2),1))/3;
    ylim([tmp1, tmp2])
    
    title(filenames{i},'interpret','none')
end


%% tSNE
figure(3)
data = pooled_data(:,1:end-1);
labels = pooled_data(:,end);
mappedX = fast_tsne(data); % run the fast bh version of tsne, downloaded from Laurens van der Maaten's website
subplot(1,2,1); plot(mappedX(:,1),mappedX(:,2),'g.')
subplot(1,2,2); FlowJo_contour2D(mappedX(:,1),mappedX(:,2),1)
figure(4)
counter = 1;
for i=[8 11 14 62 63 65]  % visualize 6 selected samples from the total of 103 samples. These 6 belong two MOAs. The first three are Aurora Kinase inhibitors, and the last three are labeled as general Kinase inhibitors.
    subplot(2,3,counter); counter = counter+1;
    plot(mappedX(:,1),mappedX(:,2),'g.',mappedX(labels==i,1),mappedX(labels==i,2),'b.');
    title(filenames{i},'interpret','none')
end




%% SPADE 
num_clusters = 200;
data = pooled_data(:,1:end-1);
labels = pooled_data(:,end);
% cluster the cells
[idx,cluster_center,center_ind]=deterministic_kmeans_for_SPADE(data,num_clusters,ones(size(data,1),1));  
idx = idx';
cluster_center = cluster_center';
center_ind = center_ind';
% build a graph for the clusters for EMD
adj = mst_from_dist_matrix(squareform(pdist(cluster_center')));
adj_tree = adj;
KNN_parameter = 5;
ns = createns(cluster_center','nsmethod','kdtree','Distance','euclidean');
[idx_tmp, dist_tmp] = knnsearch(ns,cluster_center','k',KNN_parameter+1);
KNN_edges = unique(sort([reshape(repmat(idx_tmp(:,1),1,size(idx_tmp,2)),prod(size(idx_tmp)),1),idx_tmp(:)],2),'rows');
adj(KNN_edges(:,1)+(KNN_edges(:,2)-1)*size(adj,1))=1;
adj = double((adj + adj')~=0)-eye(size(adj));
adj_graph = adj;
% compute frequency
counts = zeros(max(labels),num_clusters);
for i=1:max(labels)
    for j=1:max(idx)
        counts(i,j) = sum(labels(:)==i & idx(:)==j);
    end
end
freq = counts./repmat(sum(counts,2),1,size(counts,2));

% draw SPADE trees
[mappedX,spring,distance]=kamada_kawai_spring_layout_mex(...
    adj_tree, 1e-30, 20000, 1, ...   % adj, tolerance, max iteration, spring_constant
    [], 1, 0, 'matrix');     % progressive_opt, options.edge_length, edge_weights, edge_weight_opt;
mappedX = mappedX';
mappedX = mappedX-repmat(mean(mappedX,2),1,size(mappedX,2));
mappedX = mappedX./max(max(abs(mappedX)))*50;
node_positions = mappedX;
figure(5)
draw_SPADE_tree_annotation(adj_tree, node_positions, 7*ones(1,num_clusters), zeros(1,num_clusters), [-1,1], 1, 0, zeros(1,num_clusters), 'jet', [],[]);
title('SPADE tree')


% visualize 6 selected samples from the total of 103 samples. The same set of 6 as above
figure(6)
counter = 1;
for i=[8 11 14 62 63 65]
    subplot(2,3,counter); counter = counter+1;
    draw_SPADE_tree_annotation(adj_tree, node_positions, 7*ones(1,num_clusters), freq(i,:), max(freq(i,:))*[-1,1], 1, 0, zeros(1,num_clusters), 'jet', [],[]);
    title(filenames{i},'interpret','none')
end



%% define distance measure - corr - sample level visualization
figure(7)
% heatmap and hierarchical clustering
pairwise_dist_corr = 1-show_sample_corr(freq');
hierarchical_clustering_links = linkage(pairwise_dist_corr,'single');
[H,T,perm_r] = dendrogram(hierarchical_clustering_links,size(pairwise_dist_corr,1));
perm_c = perm_r;
h=subplot(11,2,[1:2:16]);
imagesc(pairwise_dist_corr(perm_r,perm_c));
title('Hierarchial clustering of pairwise Corr distance of SPADE derived distributions')
h=subplot(11,2,[19]);
[~,~,label_moa] = unique(MOA);
imagesc(label_moa(perm_r)')

% SPADE-like MST visualization
adj_samples = mst_from_dist_matrix(pairwise_dist_corr);
[mappedX,spring,distance]=kamada_kawai_spring_layout_mex(...
    adj_samples, 1e-30, 20000, 1, ...   % adj, tolerance, max iteration, spring_constant
    [], 1, 0, 'matrix');     % progressive_opt, options.edge_length, edge_weights, edge_weight_opt;
mappedX = mappedX';
mappedX = mappedX-repmat(mean(mappedX,2),1,size(mappedX,2));
mappedX = mappedX./max(max(abs(mappedX)))*50;
h=subplot(11,2,[2:2:16]);
scatter(mappedX(1,:),mappedX(2,:),80,label_moa,'fill')
for i=1:length(label_moa)
    text(mappedX(1,i)+0.75,mappedX(2,i),num2str(label_moa(i)))
end
title('MST based on pairwise Corr distance of SPADE derived distributions')
axis off
h=subplot(11,2,[20]);
imagesc(label_moa(:)')
set(h,'ytick',[],'xtick',[])
for i=1:12
    text(mean(find(label_moa==i)),1,num2str(i));
    text(mean(find(label_moa==i)),2,MOA{find(label_moa==i,1)},'rotation',-60);
end
colormap jet

    
%% define distance measure - EMD
[is_connected, shortest_hop] = check_graph_connectivity(adj_graph.*squareform(pdist(cluster_center')));
fprintf('Computing pairwise EMDs for %d samples: %5d,%5d',size(freq,1),0,0);
EMD = zeros(size(freq,1));
for i=1:size(freq,1)-1
    for j=i+1:size(freq,1)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d,%5d',i,j);
        EMD(i,j) = EarthMoverDist2(freq(i,:),freq(j,:),shortest_hop-diag(diag(shortest_hop))); 
        drawnow;
    end
end
pairwise_dist_EMD = EMD + EMD';


figure(8)
hierarchical_clustering_links = linkage(pairwise_dist_EMD,'single');
[H,T,perm_r] = dendrogram(hierarchical_clustering_links,size(pairwise_dist_EMD,1));
perm_c = perm_r;
h=subplot(11,2,[1:2:16]);
imagesc(pairwise_dist_EMD(perm_r,perm_c));
title('Hierarchial clustering of pairwise EMD distance of SPADE derived distributions')
h=subplot(11,2,[19]);
[~,~,label_moa] = unique(MOA);
imagesc(label_moa(perm_r)')

% SPADE-like MST visualization
adj_samples = mst_from_dist_matrix(pairwise_dist_EMD);
[mappedX,spring,distance]=kamada_kawai_spring_layout_mex(...
    adj_samples, 1e-30, 20000, 1, ...   % adj, tolerance, max iteration, spring_constant
    [], 1, 0, 'matrix');     % progressive_opt, options.edge_length, edge_weights, edge_weight_opt;
mappedX = mappedX';
mappedX = mappedX-repmat(mean(mappedX,2),1,size(mappedX,2));
mappedX = mappedX./max(max(abs(mappedX)))*50;
h=subplot(11,2,[2:2:16]);
scatter(mappedX(1,:),mappedX(2,:),80,label_moa,'fill')
for i=1:length(label_moa)
    text(mappedX(1,i)+0.75,mappedX(2,i),num2str(label_moa(i)))
end
title('MST based on pairwise EMD distance of SPADE derived distributions')
axis off
h=subplot(11,2,[20]);
imagesc(label_moa(:)')
set(h,'ytick',[],'xtick',[])
for i=1:12
    text(mean(find(label_moa==i)),1,num2str(i));
    text(mean(find(label_moa==i)),2,MOA{find(label_moa==i,1)},'rotation',-60);
end
colormap jet


%% compare corr and EMD based sample level visualization
figure(9); 
plot(pairwise_dist_corr(:),pairwise_dist_EMD(:),'.')
xlabel('pairwise corr distance')
ylabel('pairwise EMD distance')
