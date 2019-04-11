function model = vbr(X, y)
%% pre-compute some constants
maxiter = 100; tol=1e-3;
[M, N] = size(X);
lam = 0;
theta = ones(N,1);
g0 = 1e-2*ones(N,1); h0 = 1e-4*ones(N,1);
c0 = 1e-2; d0 = 1e-4;
e0 = 1; f0 = N;
g = g0 + 0.5;
c = c0 + 0.5*M;
e = e0; f = f0;
E_lambda = g0 ./ h0;
E_tau = c0 / d0;
X2 = X'*X;
Xy = X'*y;
y2 = y'*y;

[theta, mu, E_lambda,E_tau,g,e,f] = update(theta,X2,E_lambda,E_tau,Xy,h0,g,d0,y2,c,X,y,e,f,e0,f0);
MU = mu;

converge = false; t = 1; 
while ~converge
    % update posterior parameters of a based on xi
    [theta, mu_new, E_lambda,E_tau,g,e,f] = update(theta,X2,E_lambda,E_tau,Xy,h0,g,d0,y2,c,X,y,e,f,e0,f0);

    MU = [MU, mu_new];
    % converge or not
    if t>= maxiter
        converge = true;
    elseif sum(abs(mu_new - mu)) < tol
        converge = true;
    end
    %ELBO_last = ELBO(t);
    t = t + 1;
    mu = mu_new;
end
gamma = theta;
gamma(gamma<0.5) = 0; gamma(gamma>=0.5) = 1;
gamma = logical(gamma);
p1 = sum(gamma~=0);
coef = zeros(size(gamma));
coef(gamma) = mu(gamma);
BIC = 2 * sum(log(1+exp( - y .* (X * coef)))) + p1 * log(M);

%% output


model.BIC = BIC;
model.mu = mu;
model.coef = coef;
model.Step = t;
model.theta = theta;
model.lambda = lam;
model.muHistory = MU;
model.gamma = gamma;

model.converge_flag = (t < maxiter);
end

function [theta, mu, E_lambda,E_tau,g,e,f] = update(theta,X2,E_lambda,E_tau,Xy,h0,g,d0,y2,c,X,y,e,f,e0,f0)
Omega = theta*theta'+diag(theta.*(1-theta));
X2Omega = X2 .*Omega;
% b
invSigma = diag(E_lambda) + E_tau * X2Omega;
Sigma = inv(invSigma);
mu = E_tau * Sigma * bsxfun(@times, Xy, theta);
D_w = Sigma + mu * mu';
% lambda
h = h0 + 0.5 * diag(D_w);
E_lambda = g./h;
% tau
d = d0 + 0.5*(y2-2*Xy'*(mu.*theta)+trace(D_w*(X2Omega)));
E_tau = c/d;
% a
temp1 = X'*(y-X*(mu.*theta)).*mu;
temp2 = diag(X2) .* (mu.^2.*theta-0.5*diag(D_w));
temp = psi(e)-psi(f)+E_tau*(temp1+temp2);
theta = 1./(1+exp(-temp));
theta(theta>0.9999)=0.9999; theta(theta<0.0001)=0.0001;
% rho
e = e0 + sum(theta);
f = f0 + sum(1-theta);
end 