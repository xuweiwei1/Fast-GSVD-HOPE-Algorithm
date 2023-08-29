function [U,V]=fast_embed_main(A,K,beta)
% Input: 
% A: N*N adjacency matrix (sparse)
% K: dimensionality of embedding space
% beta: 4
% Output:
% U: N*K left embedding matrix
% V: N*K right embedding matrix
% The high-order proximity (katz) matrix is approximated by U * V'

[N, ~] = size(A);

% Katz: S = sum_{l=1}^{+inf}{beta*A}^l
if nargin < 3
    beta = 0.5 / getRadius(A);
end

A=full(A);
Mb = beta.*A;
Ma = eye(N)-Mb;

H=[Ma';Mb'];
[L,~]=qr(H,0);
L1=L(1:N,:);
L2=L(N+1:end,:);
[S,psi]=eig(L1'*L1);
phi=sqrt(eye(N)-psi);
psi=sqrt(psi);
u=L1*S;
v=L2*S;
U=zeros(n,r);
V=zeros(n,r);
for i=1:K
    U(:,i)=u(:,i)/psi(i,i)*sqrt(phi(i,i)/psi(i,i));
    V(:,i)=v(:,i)/phi(i,i)*sqrt(phi(i,i)/psi(i,i));
end

end