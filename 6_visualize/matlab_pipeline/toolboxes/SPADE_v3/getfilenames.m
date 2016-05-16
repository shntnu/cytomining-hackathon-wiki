function filenames = getfilenames(is_order_files)
% is_order_files can be a string of folder name
% if is_order_files is not defined, just get the filenames
% if is_order_files==1, get filenames and order by length of names
% if is_order_files~=1, just get the filenames, no ordering

if exist('is_order_files') && isstr(is_order_files)
    tmp = dir(is_order_files);
else
    tmp = dir;
end
filenames = cell(length(tmp),1);
for i=1:length(filenames)
    filenames{i} = tmp(i).name;
end

filenames(ismember(filenames,{'.'}) | ismember(filenames,{'..'}))=[];

if exist('is_order_files') && isnumeric(is_order_files) && is_order_files==1
    len = zeros(1,length(filenames));
    for i=1:length(filenames)
        len(i) = length(filenames{i});
    end
    [Y,I] = sort(len);
    filenames = filenames(I);
end
