function [Sigma]=FuncICP_crossCovar(P,X,mu_p,mu_x)
%This function calculates the cross-covariance matrix Sigma from the data
%sets P and X.

[N,M]=size(P);
Sigma_temp=zeros(N);

for i=1:M
    Sigma_temp=Sigma_temp+P(:,i)*X(:,i)'-mu_p*mu_x';
end

Sigma=Sigma_temp/M;


