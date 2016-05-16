function [numLines] = get_txt_num_lines_slow(filename)


fid = fopen(filename);
tmp = fgetl(fid); line_count = 1; 
while ischar(tmp) 
    tmp = fgetl(fid); line_count = line_count+1;
end
fclose(fid);
line_count = line_count-1;
numLines = line_count;

return


