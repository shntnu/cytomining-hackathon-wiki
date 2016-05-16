function [data] = flow_arcsinh(data,cofactor)
% cofactor can be one signal number, for all channels
%          can be a vector, one element for one channel
% if cofactor == 0, then it means linear, no transform. 

if length(cofactor)~=1 && length(cofactor)~=size(data,1)
    disp('wrong input of cofactors')
    disp('data not transformed')
    return
end

if length(cofactor)==1
    if cofactor==0
        return
    else
        data = log( data(:,:)/cofactor + sqrt((data(:,:)/cofactor).^2+1) );
        return
    end
end

if length(cofactor) == size(data,1)
    for i=1:size(data,1)
        if cofactor(i)==0
            continue;
        else
            data(i,:) = log( data(i,:)/cofactor(i) + sqrt((data(i,:)/cofactor(i)).^2+1) );
        end        
    end
end

