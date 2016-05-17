function [pairs] = find_matrix_big_element(matrix, threshold)

ind = find(matrix>=threshold);
j = ceil(ind/size(matrix,1));
i = ind - (j-1)*size(matrix,1);
pairs = [i,j];