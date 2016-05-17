function [numLines] = get_txt_num_lines_fast(filename)
    fid = fopen(filename);
    d = fread(fid,[1,inf],'char'); 
    fclose(fid);
    numLines = sum(d==10);
    if d(end)~=10
        numLines = numLines + 1;
    end

return
