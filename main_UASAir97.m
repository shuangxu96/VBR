task = 'USAir97';
load(['E:\code\network_reconstruction\power\',task,'.mat'])
A=full(Problem.A); A(A~=0)=1;
M = 100;
N = length(A);
maxiter = 10;

for numiter=1:maxiter
    display(numiter)
    % simulate ele model 
    [I, V, W] = ElectricalCurrent(A,M);
    % Reconstruct network
    T_cv = zeros(size(A)); temp_cv = zeros(1,N);
    T_vb = zeros(size(A)); temp_vb = zeros(1,N); 
    for i=1:N
        G = V - V(:,i);
        % cv+lasso
        tic;[B,info] = lasso(G,I(:,i),'CV',5);temp_cv(i)=toc;
        T_cv(:,i) = B(:,info.Index1SE);
        % bayes recon
        G(:,i)=[];
        tic;model = vbr(G,I(:,i));temp_vb(i)=toc;
        temp = model.theta';
        T_vb(:,i) = [temp(1:(i-1)),0,temp(i:end)];
        temp = model.coef';
        W_vb(:,i) = [temp(1:(i-1)),0,temp(i:end)];
    end
	
    MSE_cv(numiter)=sqrt(sum(sum ((T_cv-W).^2) )) / sqrt(sum(W(:).^2));
    MSE_vb(numiter)=sqrt(sum(sum ((W_vb-W).^2) )) / sqrt(sum(W(:).^2));
    
    T_cv2=T_cv;T_cv2(T_cv2~=0)=1;
    T_vb(T_vb>0.5)=1;T_vb(T_vb<=0.5)=0;

    time_cv(numiter) = sum(temp_cv); 
	time_vb(numiter) = sum(temp_vb); 
    [TP_cv(numiter), TN_cv(numiter)] = Metrics(A, T_cv2);
    [TP_vb(numiter), TN_vb(numiter)] = Metrics(A, T_vb);
    
    [[mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
    [mean(TP_cv),std(TP_cv),mean(TP_vb),std(TP_vb)]
    [mean(TN_cv),std(TN_cv),mean(TN_vb),std(TN_vb)]
    [mean(MSE_cv),std(MSE_cv),mean(MSE_vb),std(MSE_vb)]]
    save([task,'_ECurrent_M',num2str(M),'.mat'])
end
result1 = [[mean(TP_cv),std(TP_cv);mean(TP_vb),std(TP_vb)], ...
			[mean(TN_cv),std(TN_cv);mean(TN_vb),std(TN_vb)],...
			[mean(MSE_cv),std(MSE_cv); mean(MSE_vb),std(MSE_vb)],...
			[mean(time_cv),std(time_cv);mean(time_vb),std(time_vb)]]
			
M=400;
for numiter=1:maxiter
    display(numiter)
    % Simulate Communication model 
    [F, O, W] = Communication(A,M);
    % Reconstruct network
    T_cv = zeros(size(A)); temp_cv = zeros(1,N);
    T_vb = zeros(size(A)); temp_vb = zeros(1,N); 
    for i=1:N
        % cv+lasso
		tic;[B,info] = lasso(O,F(:,i),'CV',5);temp_cv(i)=toc;
        T_cv(:,i) = B(:,info.Index1SE);
        % bayes recon
        tic;model = vbr(O,F(:,i));temp_vb(i)=toc;
        T_vb(:,i) = model.theta;%[temp(1:(i-1)),0,temp(i:end)];
        W_vb(:,i) = model.coef;
    end
	
    MSE_cv(numiter)=sqrt(sum(sum ((T_cv-W).^2) )) / sqrt(sum(W(:).^2));
    MSE_vb(numiter)=sqrt(sum(sum ((W_vb-W).^2) )) / sqrt(sum(W(:).^2));
    
    T_cv2=T_cv;T_cv2(T_cv2~=0)=1;
    T_vb(T_vb>0.5)=1;T_vb(T_vb<=0.5)=0;

    time_cv(numiter) = sum(temp_cv); 
	time_vb(numiter) = sum(temp_vb); 
    [TP_cv(numiter), TN_cv(numiter)] = Metrics(A, T_cv2);
    [TP_vb(numiter), TN_vb(numiter)] = Metrics(A, T_vb);
    
    [[mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
    [mean(TP_cv),std(TP_cv),mean(TP_vb),std(TP_vb)]
    [mean(TN_cv),std(TN_cv),mean(TN_vb),std(TN_vb)]
    [mean(MSE_cv),std(MSE_cv),mean(MSE_vb),std(MSE_vb)]]
    save([task,'_Commnu_M',num2str(M),'.mat'])
end
result2 = [[mean(TP_cv),std(TP_cv);mean(TP_vb),std(TP_vb)], ...
			[mean(TN_cv),std(TN_cv);mean(TN_vb),std(TN_vb)],...
			[mean(MSE_cv),std(MSE_cv); mean(MSE_vb),std(MSE_vb)],...
			[mean(time_cv),std(time_cv);mean(time_vb),std(time_vb)]]
% PLplot(A)
% fig=gcf
% fig.PaperSize=[5,5]
% fig.InnerPosition=[400,400,420,380]
% saveas(gcf,'USAir','epsc')