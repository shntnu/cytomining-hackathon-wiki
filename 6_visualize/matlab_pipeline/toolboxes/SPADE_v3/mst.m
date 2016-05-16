function [adj,adj2, cost_value] = mst(X,working_mode,exclude_adj)

% Minimal or Minimum Spanning Tree based on Euclidian distances
% MST in short: use (X)(n x p) to form (n-1) lines to connect (n) objects in the shortest possible way in the (p)
% dimensional variable-space, under the condition 'no closed loops allowed'.
% working_mode : 'euclidean' (default)
%                'corr'
%                'abs_corr'
% 
% out:Xmst (objects-1 x 2) link set between 'objects' indexed as rows in X
%     adj adjacency matrix


cost_value = 0;

if ~exist('working_mode')
    working_mode = 'euclidean';
end
if isempty(intersect({'euclidean','corr','abs_corr'}, working_mode))
    working_mode = 'euclidean';
end
if ~isempty(intersect({'corr','abs_corr'}, working_mode))
    X = per_gene_normalization(X); % this makes computing correlation easier
    X = X./norm(X(1,:));
end

if ~exist('exclude_adj') || isempty(exclude_adj) 
    exclude_adj = sparse(size(X,1),size(X,1));
end

[nX,mX] = size(X);

components = []; active_components =[];
adj = sparse(size(X,1),size(X,1)); 
adj2 = sparse(size(X,1),size(X,1)); 
count = 0; %fprintf('constructing a total of %d MST edges ... %6d', size(X,1)-1,count);
for i=1:nX
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(X,i,1:nX); 
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(X,i,1:nX);
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(X,i,1:nX); 
    end
    dist(i) = max(dist)+1; 
    dist = dist + exclude_adj(i,:).*(max(dist)+1);    
    dist = full(dist);
    [Dmin,Dwin] = min(dist);
    Xmst(i,:) = [i Dwin];
    if adj(i,Dwin)==0 && adj(Dwin,i)==0
        adj(i,Dwin)=1;adj(Dwin,i)=1;
        adj2(i,Dwin)=Dmin;adj2(Dwin,i)=Dmin;
        cost_value = cost_value + Dmin;
        count = count + 1; %fprintf('\b\b\b\b\b\b%6d', count);
    end
    if isempty(components)
        components = sparse(zeros(size(X,1),1)); components([i,Dwin],1) = 1; active_components=1;
    else
        [existing_comp1] = find(components(i,:)==1 & active_components==1);
        [existing_comp2] = find(components(Dwin,:)==1 & active_components==1);
        if isempty(existing_comp1) && isempty(existing_comp2)
            components = [components,zeros(size(X,1),1)]; components([i,Dwin],end) = 1; active_components = [active_components,1];
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
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(X,ind1,ind2); 
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(X,ind1,ind2);
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(X,ind1,ind2);
    end
    dist = dist + exclude_adj(ind1,ind2).*(max(max(dist))+1); dist = full(dist);
    [Dmin,ind] = min(reshape(dist,length(ind1)*length(ind2),1));
    j = ceil(ind/length(ind1));
    i = ind - (j-1)*length(ind1);
    Xmst = [Xmst; [ind1(i),ind2(j)]];
    adj(ind1(i),ind2(j))=1; adj(ind2(j),ind1(i))=1; 
    adj2(ind1(i),ind2(j))=Dmin; adj2(ind2(j),ind1(i))=Dmin; 
    cost_value = cost_value + Dmin;
    [existing_comp2] = find(components(ind2(j),:)==1 & active_components==1);
    components(:,existing_comp1) = components(:,existing_comp1) + components(:,existing_comp2);
    active_components(existing_comp2)=0;
    count = count + 1; %fprintf('\b\b\b\b\b\b%6d', count);
end
% fprintf('\n');
return




function dist = comp_dist_euclidean(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
for i=1:length(ind1)
    dist(i,:) = sqrt(sum((repmat(X(ind1(i),:),length(ind2),1) - X(ind2,:)).^2,2)); 
end
return


function dist = comp_dist_corr(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
corr = X(ind1,:)*X(ind2,:)';
dist = 1-corr; 
return


function dist = comp_dist_abs_corr(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
corr = X(ind1,:)*X(ind2,:)';
dist = 1-abs(corr); 
return

