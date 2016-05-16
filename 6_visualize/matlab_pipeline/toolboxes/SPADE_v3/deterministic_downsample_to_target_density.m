function [ind_selected] = deterministic_downsample_to_target_density(data, local_density, target_density)

keep_prob = min(1,(target_density./local_density));

N = length(local_density);
ind_available = 1:N;
ind_selected = zeros(1,N);

fprintf('deterministic downsampling in progress ... %3d%%',0);
while ~isempty(ind_available)
    [~,i] = max(local_density);
    ind_selected(i)=1;
    dist = sum(abs(bsxfun(@minus,data(:,i),data(:,ind_available))),1);
    [~,I] = sort(dist,'ascend');
    [k] = max(sum(cumsum(keep_prob(ind_available(I)))<=1),1);
    keep_prob(ind_available(I(1:k)))=0;
    local_density(ind_available(I(1:k)))=0;
    ind_available(I(1:k))=[];
    fprintf('\b\b\b\b%3d%%',floor((1-length(ind_available)/N)*100));
%     figure(2)
%     plot(data(1,:),data(2,:),'.b',data(1,setdiff(1:N,ind_available)),data(2,setdiff(1:N,ind_available)),'.g',data(1,ind_selected==1),data(2,ind_selected==1),'o');
%     drawnow;
end




% keep_prob = min(1,(target_density./local_density));
% 
% ind_available = ones(1,length(local_density));
% ind_selected = zeros(1,length(local_density));
% 
% fprintf('deterministic downsampling in progress ... %3d%%',0);
% while max(ind_available)~=0
%     [~,i] = max(local_density);
%     ind_selected(i)=1;
%     dist = sum(abs(bsxfun(@minus,data(:,i),data)),1);
%     [~,I] = sort(dist,'ascend');
%     [k] = max(sum(cumsum(keep_prob(I))<=1),1);
%     keep_prob(I(1:k))=0;
%     ind_available(I(1:k))=0;
%     local_density(I(1:k))=0;
%     fprintf('\b\b\b\b%3d%%',floor((1-mean(ind_available))*100));
% %     plot(data(1,:),data(2,:),'.b',data(1,ind_available==0),data(2,ind_available==0),'.g',data(1,ind_selected==1),data(2,ind_selected==1),'o');
% %     drawnow;
% end
% 







