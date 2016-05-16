function table_content = read_xls_v2(filename, sep, skip_rows, column_ind, num_rows)
% table_content = read_xls_v2(filename, sep, skip_rows, column_ind, num_rows)

table_content=[];
if exist('num_rows') 
    total_lines = skip_rows + num_rows;
else
    try
        total_lines = get_txt_num_lines_fast(filename);
    catch
        total_lines = get_txt_num_lines_slow(filename);
    end
end

if ~exist('column_ind')
    column_ind = [];
end

fid = fopen(filename);
ftext = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', skip_rows+1);
if iscell(ftext{1})
    tmp = ftext{1}{1}; 
else
    tmp = ftext{1};
end
while tmp(1)==sep, tmp = tmp(2:end); end % the first char can not be sep
if tmp(end)~=sep, tmp = [tmp, sep]; end % the last char has to be sep    
num_columns = sum( tmp==sep);
fclose(fid);

chunkSize = 5000;
if total_lines-skip_rows <=chunkSize
    fid = fopen(filename);
    if skip_rows~=0
        ftext = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', skip_rows-1); % skip the "skip_rows"
    end
    ftext = textscan(fid, '%s', (total_lines-skip_rows)*num_columns, 'delimiter', sep);
    if isempty(column_ind)
        table_content = reshape(ftext{1},num_columns,total_lines-skip_rows)';
    else
        table_content = reshape(ftext{1},num_columns,total_lines-skip_rows)';
        table_content = table_content(:,column_ind);
    end
    fclose(fid);
else
    fid = fopen(filename);
    if skip_rows~=0
        ftext = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', skip_rows-1); % skip the "skip_rows"
    end
    lines_left_to_read = total_lines - skip_rows;
    while 1
        if lines_left_to_read<=0  % all needed lines are read already
            break;
        end
        this_chunkSize = min(chunkSize, lines_left_to_read);
        ftext = textscan(fid, '%s', this_chunkSize*num_columns, 'delimiter', sep);
%         tmp = parse_txt_lines(ftext{1},sep,column_ind);
        if isempty(column_ind)
            table_content = [table_content; reshape(ftext{1},num_columns,this_chunkSize)'];
        else
            tmp =  reshape(ftext{1},num_columns,this_chunkSize)';
            tmp = tmp(:,column_ind);
            table_content = [table_content;tmp];
        end
        lines_left_to_read = lines_left_to_read - this_chunkSize
    end
    fclose(fid);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table_content = parse_txt_lines(texts,sep,column_ind)

% get the total number of seq points in the first 10 lines
len_separation_points = 0;
for i=1:min(10,length(texts))
    tmp = texts{i};
    if tmp(end)~=sep, tmp=[tmp,sep]; end
    separation_points = [0,find(tmp==sep)];
    len_separation_points = max(len_separation_points, length(separation_points));
end
if isempty(column_ind)
    column_ind = 1:(len_separation_points-1);
end


table_line_count = 1;
for i=1:length(texts)
    tmp = texts{i};
    if length(tmp)>=24 && sum(tmp(1:24)=='!series_matrix_table_end')==24 % GEO serial matrix last line usually look like this
        break;
    end
    while isempty(tmp) || length(tmp)<length(separation_points)-1 % if this line is empty, skip it
        continue;
    end
    if tmp(end)~=sep, tmp=[tmp,sep]; end
    separation_points = [0,find(tmp==sep)];

    % now focus on this particular line and parse this line
    for j=1:length(separation_points)-1
        if sum(column_ind==j)==1
            seg = tmp(separation_points(j)+1:separation_points(j+1)-1);
            col = find(column_ind==j);
%             if isanumber(seg)==1 
%                 table_content(table_line_count,col) = {str2num(seg)};
%             else
                table_content(table_line_count,col) = {seg};
%             end
        end
    end
    table_line_count = table_line_count + 1;
end

return






    
function [result] = isanumber(str)
if isempty(str)
    result=0; return
end
while (str(1)==char(9) || str(1)==32 || str(1)==13 || str(1)==',') && length(str)>=2 % 9--tab, 32-space, 13-enter
    str = str(2:end);
end
while (str(end)==char(9) || str(end)==32 || str(end)==13 || str(end)==',') && length(str)>=2 % 9--tab, 32-space, 13-enter
    str = str(1:end-1);
end
if ~isempty(setdiff(str, '0123456789.+-E'))
    result = 0; return
end
if ~isempty(intersect(str,'+-')) && ~isempty(setdiff(find(str=='+' | str=='-'),1)) && ~isempty(setdiff(str(setdiff(find(str=='+' | str=='-'),1)-1),'eE'))
    result=0; return
end
if length(str2num(str))~=1
    result = 0; return
end
result = 1;
return
