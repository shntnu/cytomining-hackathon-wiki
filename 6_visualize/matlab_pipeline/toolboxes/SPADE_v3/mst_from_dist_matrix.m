function [adj,adj2, cost_value] = mst_from_dist_matrix(dist_matrix)
% [adj,adj2, cost_value] = mst(dist_matrix)
% Minimal or Minimum Spanning Tree based on dist_matrix

cost_value = 0;
[nX] = size(dist_matrix,1);
dist_matrix = dist_matrix+eye(size(dist_matrix))*max(max(dist_matrix))*2; % make sure the diag element are high, so that we don't have self-loops

components = []; active_components =[];
adj = sparse(nX,nX); 
adj2 = sparse(nX,nX); 

fprintf('creating a total of %d edges ... %6d',nX-1,0); count = 0;
for i=1:nX
    dist = dist_matrix(i,:);
    [Dmin,Dwin] = min(dist);
    Xmst(i,:) = [i Dwin];
    if adj(i,Dwin)==0 && adj(Dwin,i)==0
        adj(i,Dwin)=1;adj(Dwin,i)=1;
        if Dmin==0, Dmin = 1e-10; end
        adj2(i,Dwin)=Dmin;adj2(Dwin,i)=Dmin;
        cost_value = cost_value + Dmin;
        count = count +1; fprintf('\b\b\b\b\b\b%6d',count);
    end
    if isempty(components)
        components = sparse(zeros(nX,1)); components([i,Dwin],1) = 1; active_components=1;
    else
        [existing_comp1] = find(components(i,:)==1 & active_components==1);
        [existing_comp2] = find(components(Dwin,:)==1 & active_components==1);
        if isempty(existing_comp1) && isempty(existing_comp2)
            components = [components,zeros(nX,1)]; components([i,Dwin],end) = 1; active_components = [active_components,1];
        elseif ~isempty(existing_comp1) && isempty(existing_comp2)
            components([i,Dwin],existing_comp1)=1;
        elseif isempty(existing_comp1) && ~isempty(existing_comp2)
            components([i,Dwin],existing_comp2)=1;
        elseif ~isempty(existing_comp1) && ~isempty(existing_comp2) && existing_comp1~=existing_comp2
            components = [components, components(:,existing_comp1)+components(:,existing_comp2)];
            active_components = [active_components,1];
            active_components([existing_comp1, existing_comp2])=0;
        end
    end
end
 
while sum(active_components)>1
%     sum(active_components)
    components_sizes = sum(components); components_sizes(active_components==0) = max(components_sizes+1);
    [dummy, existing_comp1] = min(components_sizes);
    ind1 = find(components(:,existing_comp1)==1); ind1 = ind1(:)';
    ind2 = setdiff(1:size(components,1),ind1); ind2 = ind2(:)';
    dist = dist_matrix(ind1,ind2);
    [Dmin,ind] = min(reshape(dist,length(ind1)*length(ind2),1));
    j = ceil(ind/length(ind1));
    i = ind - (j-1)*length(ind1);
    Xmst = [Xmst; [ind1(i),ind2(j)]];
    adj(ind1(i),ind2(j))=1; adj(ind2(j),ind1(i))=1; 
    if Dmin==0, Dmin = 1e-10; end
    adj2(ind1(i),ind2(j))=Dmin; adj2(ind2(j),ind1(i))=Dmin; 
    cost_value = cost_value + Dmin;
    [existing_comp2] = find(components(ind2(j),:)==1 & active_components==1);
    components(:,existing_comp1) = components(:,existing_comp1) + components(:,existing_comp2);
    active_components(existing_comp2)=0;
    count = count +1; fprintf('\b\b\b\b\b\b%6d',count);
end
fprintf('\n');

return


