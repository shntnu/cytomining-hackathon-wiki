function [perm_r,perm_c] = HCC_heatmap(data,display_mode)


% This part is for visualization the training data using hierarchical clustering provided by the matlab statistical toolbox
if exist('display_mode')
    if display_mode=='r'
        eucD = pdist(data,'euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_r] = dendrogram(clustTreeEuc,size(data,1));
        perm_c = 1:size(data,2);
    elseif display_mode=='c'
        eucD = pdist(data','euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2));
        perm_r = 1:size(data,1);
    else
        eucD = pdist(data,'euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_r] = dendrogram(clustTreeEuc,size(data,1));
        eucD = pdist(data','euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2));
    end
    imagesc(data(perm_r,perm_c));
    return
end

eucD = pdist(data,'euclidean');
clustTreeEuc = linkage(eucD,'average');
[H,T,perm_r] = dendrogram(clustTreeEuc,size(data,1));
eucD = pdist(data','euclidean');
clustTreeEuc = linkage(eucD,'average');
[H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2));
imagesc(data(perm_r,perm_c));
