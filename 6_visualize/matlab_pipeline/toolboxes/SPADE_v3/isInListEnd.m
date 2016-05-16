function ind = isInListEnd(all_filenames,endstring)
% function ind = isInList(all_filenames,endstring)
% 

all_filenames = upper(all_filenames);
endstring = upper(endstring);
ind = [];
for i=1:length(all_filenames)
    if length(all_filenames{i})>=length(endstring) && isequal(all_filenames{i}(end+1-length(endstring):end),endstring)
       ind = [ind,i];
    end
end

