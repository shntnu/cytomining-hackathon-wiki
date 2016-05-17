function [corr_matrix] = show_sample_corr(data)

if sum(sum(isnan(data)))==0
    data = per_gene_normalization(data')';
    data = data/norm(data(:,1));
    corr_matrix = data'*data;
else
    corr_matrix = eye(size(data,2),size(data,2));
    for i=1:size(data,2)-1
        for j=i+1:size(data,2)
            tmp = data(:,[i,j]);
            tmp(sum(isnan(tmp),2)~=0,:)=[];
            if size(tmp,1)<=2
                corr_matrix(i,j)=0;
                corr_matrix(j,i)=0;
            else
                tmp = per_gene_normalization(tmp')';
                tmp = tmp/norm(tmp(:,1));
                corr_matrix([i,j],[i,j]) = tmp'*tmp;
            end
        end
    end
end
