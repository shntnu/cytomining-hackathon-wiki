function check_marker_name_consistency

filenames = getfilenames;
filenames = filenames(isInListEnd(upper(filenames),upper('.fcs')));

all_marker_names = [];
all_channel_names = [];
for i=1:length(filenames)
    [data, marker_names, channel_names] = readfcs(filenames{i});
    all_marker_names = [all_marker_names; marker_names(:)'];
    all_channel_names = [all_channel_names; channel_names(:)'];
end

problem_with_marker_names = 0;
for i=1:size(all_marker_names,2)
    if length(unique(all_marker_names(:,i)))~=1
        problem_with_marker_names  = 1;
    end
end

problem_with_channel_names = 0;
for i=1:size(all_channel_names,2)
    if length(unique(all_channel_names(:,i)))~=1
        problem_with_channel_names  = 1;
    end
end

if problem_with_marker_names==0 && problem_with_channel_names==0
    display('all marker channel names are consistent. Great!')
    return
end


header = [{'filen\marker'}];
for i=1:size(all_marker_names,2), header = [header,{i}]; end
left_column = [];
for i=1:size(all_marker_names,1), left_column = [left_column;{i}]; end

[left_column, filenames(:)]
[header; left_column, all_marker_names]
[header; left_column, all_channel_names]
