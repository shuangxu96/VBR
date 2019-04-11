close all
clear
clc

N = 50;
M = 150;
nettype = 'ws'; % BA or WS network 

NumNoiseScale = 10;
result        = cell(1,NumNoiseScale);
maxiter       = 100;
TP_cv = zeros(1,maxiter); time_cv = TP_cv; TN_cv = TP_cv; MSE_cv = TP_cv; 
TP_vb = zeros(1,maxiter); time_vb = TP_vb; TN_vb = TP_vb; MSE_vb = TP_vb; 

for kk=1:NumNoiseScale
    scale=kk/10;
for numiter=1:maxiter
    display(['Number of Iteration:',num2str(numiter)])
    display(['Nodes:',num2str(N),'Observations',num2str(M),'NoiseScale',num2str(scale)])
	
	% Simulate network
    if strcmp(nettype,'ba')
        A = SFNG(N, 2);
    elseif strcmp(nettype,'ws')
        A = full(adjacency(WattsStrogatz(N,2,0.3)));
    end
	
    % simulate electrical current transportation model 
    [I, V, W] = ElectricalCurrent(A,M);
    I = I + randn(size(I))*scale;
	
    % Reconstruct network
    T_cv = zeros(size(A)); temp_cv = zeros(1,N);
    T_vb = zeros(size(A)); temp_vb = zeros(1,N); 
    for i=1:N
        G = V - V(:,i);
        % cv+lasso
        tic;[B,info] = lasso(G,I(:,i),'CV',5);temp_cv(i)=toc;
        T_cv(:,i) = B(:,info.Index1SE);
        % vbr
        G(:,i)=[];
        tic;model = vbr(G,I(:,i));temp_vb(i)=toc;
        temp = model.theta;
        T_vb(:,i) = [temp(1:(i-1)),0,temp(i:end)];
        temp = model.coef;
        W_vb(:,i) = [temp(1:(i-1)),0,temp(i:end)];
    end
	
	% Metrics
    MSE_cv(numiter)=sqrt(sum(sum ((T_cv-W).^2) )) / sqrt(sum(W(:).^2));
    MSE_vb(numiter)=sqrt(sum(sum ((W_vb-W).^2) )) / sqrt(sum(W(:).^2));
    
    T_cv2=T_cv;T_cv2(T_cv2~=0)=1;
    T_vb(T_vb>0.5)=1;T_vb(T_vb<=0.5)=0;

    time_cv(numiter) = sum(temp_cv); 
	time_vb(numiter) = sum(temp_vb); 
    [TP_cv(numiter), TN_cv(numiter)] = Metrics(A, T_cv2);
    [TP_vb(numiter), TN_vb(numiter)] = Metrics(A, T_vb);
    
    display([[mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
    [mean(TP_cv),std(TP_cv),mean(TP_vb),std(TP_vb)]
    [mean(TN_cv),std(TN_cv),mean(TN_vb),std(TN_vb)]
    [mean(MSE_cv),std(MSE_cv),mean(MSE_vb),std(MSE_vb)]])
end
result{kk}=[[mean(TP_cv),std(TP_cv);mean(TP_vb),std(TP_vb)], ...
			[mean(TN_cv),std(TN_cv);mean(TN_vb),std(TN_vb)],...
			[mean(MSE_cv),std(MSE_cv); mean(MSE_vb),std(MSE_vb)],...
			[mean(time_cv),std(time_cv);mean(time_vb),std(time_vb)]];
save([nettype,'_ECurrent_N',num2str(N),'_M',num2str(M),'_m',num2str(m),'.mat'],'result')
end

r=cell2mat(result');




