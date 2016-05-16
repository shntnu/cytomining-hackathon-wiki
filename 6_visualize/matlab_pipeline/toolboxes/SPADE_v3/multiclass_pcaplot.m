function [coeff,e_vectors,mean_vector,energy] = multiclass_pcaplot(data,labels,dims,plotlabels)
% [coeff,e_vectors,mean_vector,energy] = multiclass_pcaplot(data,labels,dims,plotlabels)
% each column is a point

if ~exist('dims') || isempty(dims) || min(dims)<1 || max(dims)>size(data,2)
    dims = [1,2];
end
if ~exist('plotlabels')
    plotlabels = 'o^d*p+x<>';
end


mean_vector = mean(data,2);
data = data - repmat(mean_vector,1,size(data,2));
if size(data,1)>size(data,2)
    [V,D] = eig(data'*data);
    [Y,I] = sort(diag(D),'descend');
    e_vectors = data*V(:,I(dims));
    e_vectors(:,1) = e_vectors(:,1)/norm(e_vectors(:,1));
    e_vectors(:,2) = e_vectors(:,2)/norm(e_vectors(:,2));
else
    [V,D] = eig(data*data');
    [Y,I] = sort(diag(D),'descend');
    e_vectors = V(:,I(dims));
end

coeff = pinv(e_vectors)*data;
energy = diag(D);


if ~exist('labels') | isempty(labels)
    plot(coeff(1,:),coeff(2,:),'.')
elseif ~isempty(labels)
    if ~isnumeric(labels)
        [~,~,labels] = unique(labels);
    end
    if ~isempty('plotlabels'), plotlabels = 'o^d*p+x<>'; end
    possible_labels = unique(labels);
%     for i=1:size(data,2)
%         cl = mod(find(possible_labels==labels(i)),length(plotlabels));
%         if cl==0, cl = length(plotlabels); end
%         plot(coeff(1,i),coeff(2,i),plotlabels(cl));
%         if i==1, hold on; end
%     end
%     hold off;
    hold off;
    for i=1:length(possible_labels)
        cl = mod(i,length(plotlabels));
        if cl==0, cl = length(plotlabels); end
        ind = ismember(labels,possible_labels(i));
        plot(coeff(1,ind),coeff(2,ind),plotlabels(cl));
        if i==1, hold on; end
    end
    hold off;
end

return





