function [numLines] = get_txt_num_lines(filename)


chunkSize = 100000;
numLines = 0;
fid = fopen(filename);
while ~feof(fid)
    ftext = textscan(fid, '%s', chunkSize, 'delimiter', '\n');
    numLines = numLines + numel(ftext{:});
end
fclose(fid);

return
