function [Sigma]=FuncICP_crossCovar(P,X,mu_p,mu_x)
%This function calculates the cross-covariance matrix Sigma from the data
%sets P and X.
M = size(P,2);

P_ = P - repmat(mu_p, 1, M);
X_ = X - repmat(mu_x, 1, M);
Sigma = P_*transpose(X_);


