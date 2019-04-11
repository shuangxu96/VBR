function [F,O,AW] = communication(A, M)
N = length(A(:,1));
A = double(A);
W = rand(N); 
W = bsxfun(@rdivide, W, sum(A.*W));
AW = A.*W;
O = rand(M,N)*20;
F = O*AW;
