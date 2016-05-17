function data = cell2num(table_content)

table_content_size = size(table_content);

table_content = table_content(:);
chunksize = 10000;
start_entry = 1;

fprintf('I''m working on it ... %3d%%',0);
data = zeros(size(table_content));
while start_entry<=length(data)
    end_entry = min(start_entry + chunksize - 1,  length(data));
    thisChunkSize = end_entry - start_entry + 1;
    
    tmp = str2num(cell2mat(reshape([table_content(start_entry:end_entry)'; repmat({' '},1,thisChunkSize)],1,thisChunkSize*2)));
    if length(tmp)==length(start_entry:end_entry)
        data(start_entry:end_entry) = tmp; %str2num(cell2mat(reshape([table_content(start_entry:end_entry)'; repmat({' '},1,thisChunkSize)],1,thisChunkSize*2)));
    else % there might be nulls or invalid entries
        data(start_entry:end_entry) = get_cell_number_columns(table_content(start_entry:end_entry),1,0);
    end
    start_entry = end_entry + 1;
    fprintf('\b\b\b\b%3d%%',round(end_entry/length(data)*100));
end
data = reshape(data, table_content_size);
fprintf('\b\b\b\bDone\n');
