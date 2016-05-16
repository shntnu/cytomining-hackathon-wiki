function ind = isInList(gene_names_list,gene_name)
% function ind = isInList(gene_names_list,gene_name)
% 
for i=1:length(gene_names_list)
    if isempty(gene_names_list{i})
        gene_names_list{i} = '---';
        continue;
    end
    if isnumeric(gene_names_list{i})
        gene_names_list{i} = num2str( gene_names_list{i});
    end
end

gene_names_list = upper(gene_names_list);
gene_name = upper(gene_name);
ind = [];
% ind = find(ismember(gene_names_list,gene_name)==1);
if isempty(ind)
    s= strfind(gene_names_list,gene_name);
    for i=1:length(s)
        if ~isempty(s{i})
            ind = [ind,i];
        end
    end
end

% tmp=[];
% if length(ind>1)
%     for i=1:length(ind)
%         tmp(i) = length(gene_names_list{ind(i)});
%     end
%     a = find(tmp==min(tmp));
%     a = a(1);
%     ind = ind(a);
% end
