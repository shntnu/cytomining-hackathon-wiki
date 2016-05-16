function [dist,F] = EarthMoverDist2(P,Q,C)
% dist = EarthMoverDist(P,Q,c)
% 
% distance between two probability mass distribution
% P and Q are both vectors
% Cij is the cost of moving mass between support of Pi and support of Qj
% this is to solve the following linear programming problem
%       minimize \sum_{ij} Cij*fij
%             st. fij>=0            % all earth are moving from Pi->Qj
%                 \sum_j fij = Pi   % all earth at Pi are moved away
%                 \sum_i fij = Qj   % all earch moved to Qj equals to Qj


if size(C,1)~=length(P) || size(C,2)~=length(Q)
    display('Houston, we have a problem! The dimension of inputs, P, Q, C do not match!')
    return
end


F = zeros(size(C));

% transfor the variables to the standard setting of linprog
x = F(:);
f = C(:);
A=[];b=[]; %A = -eye(length(x)); b = zeros(size(x));  % Ax<=b 
Aeq = zeros(length(P)+length(Q),prod(size(C)));
count = 1;
for i=1:length(P)
    tmp = zeros(size(C)); tmp(i,:)=1;
    Aeq(count,:) = tmp(:)'; count = count + 1;
end
for j=1:length(Q)
    tmp = zeros(size(C)); tmp(:,j)=1;
    Aeq(count,:) = tmp(:)'; count = count + 1;
end
beq = [P(:);Q(:)];
LB = zeros(size(x));
UB = ones(size(x));

options = optimset('display','off');
try
    x = linprog(f,A,b,Aeq,beq,LB,UB,'options',options);
catch
    x = linprog(f,A,b,Aeq,beq,LB,UB,[],options);
end

F = reshape(x,size(F,1),size(F,2));


dist = sum(sum(F.*C));

