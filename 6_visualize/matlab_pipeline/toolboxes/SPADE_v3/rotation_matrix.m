function R = rotation_matrix(d,v)
% R = rotation_matrix(d,v)
%     if d and v are of equal lenght, R is a rotation matrix that rotates d into v,  v = R * d;
%     otherwise, R is a rotation matrix multiplied by scaling of norm_v/norm_d, and it still satisfies v = R * d; 

initial_d = d/norm(d); % keep a back up of d
R = eye(length(d));

norm_d = norm(d);
norm_v = norm(v);

d = d/norm(d);
v = v/norm(v);

[~,ind] = sort(abs(v),'ascend');
for i=1:length(ind)-1
%     x = find((d.^2 + d(ind(i))^2).*(1-ismember((1:length(d))',ind(1:i)))>=v(ind(i)).^2,1); % wrong 
%     x = max(find((d.^2 + d(ind(i))^2).*(1-ismember((1:length(d))',ind(1:i)))>=v(ind(i)).^2));% wrong 

    x = max(setdiff(find((d.^2 + d(ind(i))^2)>=v(ind(i)).^2), ind(1:i))); % correct 
%     x = setdiff(find((d.^2 + d(ind(i))^2)>=v(ind(i)).^2), ind(1:i)); x = x(randsample(1:length(x),1)); % correct 
    y = ind(i);
    
%     [cos(t)  -sin(t)] [dx]    [sqrt(dx^2+dy^2 - vy^2]
%     [               ] [  ] =  [                     ]
%     [sin(t)   cos(t)] [dy]    [vy                   ]
    
    t = asin(  (d(x)*v(y) - d(y)*sqrt(d(x)^2+d(y)^2-v(y)^2))/(d(x)^2+d(y)^2)  );
    
    R1 = eye(length(d));
    R1(x,x) = cos(t);
    R1(y,y) = cos(t);
    R1(x,y) = sin(t);
    R1(y,x) = -sin(t);
    
    R2 = eye(length(d));
    R2(x,x) = cos(t+pi);
    R2(y,y) = cos(t+pi);
    R2(x,y) = sin(t+pi);
    R2(y,x) = -sin(t+pi);

    R3 = eye(length(d));
    R3(x,x) = cos(t);
    R3(y,y) = cos(t);
    R3(x,y) = -sin(t);
    R3(y,x) = sin(t);
    
    R4 = eye(length(d));
    R4(x,x) = cos(t+pi);
    R4(y,y) = cos(t+pi);
    R4(x,y) = -sin(t+pi);
    R4(y,x) = sin(t+pi);
    
    tmp = [abs(v - R1*d), abs(v - R2*d), abs(v - R3*d), abs(v - R4*d)];
    [~,selection] = min(tmp(y,:));
    switch selection,
        case 1 
            R = R1*R;
            d = R1*d;
        case 2
            R = R2*R;
            d = R2*d;
        case 3
            R = R3*R;
            d = R3*d;
        case 4
            R = R4*R;
            d = R4*d;
        otherwise
            1; % do nothing, this is also impossible
    end
    
%     [d v]

end

y = ind(end);
if d(y)*v(y)<0
    R1 = eye(size(R));
    R1(y,y)=-1;
    R = R1*R;
end

R = R.*(norm_v/norm_d);