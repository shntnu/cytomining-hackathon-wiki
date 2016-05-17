function [data] = get_cell_number_columns(table_content, column_ind, have_header)
% [data] = get_cell_number_columns(table_content, column_ind, have_header)

if have_header==1
    table_content = table_content(2:end,:);
end

table_content = table_content(:,column_ind);
table_config = size(table_content);

table_content = reshape(table_content, 1, prod(table_config));
data =zeros(size(table_content));

for i=1:length(data)
%     i
%     if i==998
%         i
%     end
    if isempty(table_content{i}) || (~isnumeric(table_content{i}) && isempty(str2num(table_content{i})))
        data(i) = NaN; continue;
    end
    if isnumeric(table_content{i})
        data(i) = table_content{i};
    else
        data(i) = str2num(table_content{i});
    end
end
data = reshape(data,table_config);
return
