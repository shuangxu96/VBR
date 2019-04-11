function [TPR, TNR] = Metrics(A, Ahat)
A = double(A);
P = nnz(A);
N = length(A)^2-P;
TPR = A(:)'*Ahat(:)/P;
TNR = (1-A(:))'*(1-Ahat(:))/N;

