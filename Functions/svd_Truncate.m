function [U, S, V]=svd_Truncate(T,D,th)
[U, S, V]=svd(T, 'econ');
V=V';
chi = 0;
szS = size(S,1);
for k = 1:szS
    if S(k,k)/S(1,1)>th
        chi = chi+1;
    end
end
chi = min(chi,D);
chi = max(chi,2); % avoid chi = 1
%Truncation

U = U(:,1:chi);
S = S(1:chi,1:chi);
V = V(1:chi,:);
% fprintf('The originals rank is %u and truncated rank is %u\n',szS,chi)
end