function [I, V, W] = ElectricalCurrent(A,T,param)
N = length(A(:,1));
if nargin<=2
    v0 = 1; omega = 1e+3; r = 2+rand(N);
else
    v0 = param.v0; omega = param.omega; r = param.r;
end

A = double(A).*r;
A = (A+A')/2;

V = zeros(T,N);
I = zeros(T,N);
for t=1:T
    V(t,:) = v0 * sin((omega+rand(1,N)*20)*t);
    for i=1:N
        I(t,i) = (V(t,:)-V(t,i))*A(:,i);
    end
end
W = A;
end