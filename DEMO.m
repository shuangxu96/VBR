N = 100;
M = 200;
maxiter = 100;
param.v0 = 220;
param.omega = 10^3;
param.r = rand(N)*10;
scale = 0;
m = 3;
auc_cv = zeros(1,maxiter); time_cv = auc_cv;
auc_vb = zeros(1,maxiter); time_vb = auc_vb;
for numiter=1:maxiter
    display(numiter)
    A = SFNG(N, m);
    % simulate ele model 
    [I, V] = ElectricalCurrent(A,M,param);
    I = I + randn(size(I))*scale;
    % Reconstruct network
    T_cv = zeros(size(A)); temp_cv = zeros(1,N);
    T_vb = zeros(size(A)); temp_vb = zeros(1,N); 
    for i=1:N
        G = V - V(:,i);
        % cv+lasso
        tic;[B,info] = lasso(G,I(:,i),'CV',10);temp_cv(i)=toc;
        T_cv(i,:) = B(:,info.Index1SE);
        % bayes recon
        G(:,i)=[];
        tic;model = vbr(G,I(:,i));temp_vb(i)=toc;
        temp = model.theta';
        T_vb(i,:) = [temp(1:(i-1)),0,temp(i:end)];
    end
    T_cv2=T_cv;T_cv2(T_cv2~=0)=1;
    T_vb(T_vb>0.5)=1;T_vb(T_vb<=0.5)=0;

    time_cv(numiter) = sum(temp_cv); auc_cv(numiter)=sum(sum((A-T_cv2).^2))/sum(A(:));
    time_vb(numiter) = sum(temp_vb); auc_vb(numiter)=sum(sum((A-T_vb).^2))/sum(A(:));
    [mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
    [mean(auc_cv),std(auc_cv),mean(auc_vb),std(auc_vb)]
    [median(auc_cv),std(auc_cv),median(auc_vb),std(auc_vb)]
    save(['N',num2str(N),'_M',num2str(M),'_m',num2str(m),'#',num2str(scale),'.mat'])
end
