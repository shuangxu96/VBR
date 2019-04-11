close all
clear
clc

N = 100; % Number of nodes
M = 100;
Dyn = 'C';
Alpha = -linspace(2,3,6);   % Alpha of the scale-free graph
maxiter = 100;
param.v0 = 1;
param.omega = 10^3;
param.r = 2+rand(N);

TP_cv = zeros(1,maxiter); time_cv = TP_cv; TN_cv = TP_cv; DegCorr_cv = TP_cv; 
TP_vb = zeros(1,maxiter); time_vb = TP_vb; TN_vb = TP_vb; DegCorr_vb = TP_vb; 
result = cell(1,10);
scale=0.1;
for kk=1:6
    exp_pow = Alpha(kk);
for numiter=1:maxiter
    display(uint16([N,M,scale,exp_pow,numiter]))

%define node degree distribution:
XAxis  = unique(round(logspace(0,log10(N),1000)));
YAxis  = unique(round(logspace(0,log10(N),1000))).^(exp_pow);
% create the graph with the required node degree distribution:
Graph = mexGraphCreateRandomGraph(N,XAxis,YAxis,0);
g=graph(Graph.Data(:,1),Graph.Data(:,2));
A=full(adjacency(g));
l=size(A,1);
if l<N
    A(1,(l+1):N)=1;
    A((l+1):N,1)=1;
end
    % simulate ele model 
    [I, V, W] = ElectricalCurrent(A,M,param);
    I = I + randn(size(I))*scale;
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

    time_cv(numiter) = sum(temp_cv); time_vb(numiter) = sum(temp_vb); 
    [TP_cv(numiter), TN_cv(numiter), DegCorr_cv(numiter)] = Metrics(A, T_cv2);
    [TP_vb(numiter), TN_vb(numiter), DegCorr_vb(numiter)] = Metrics(A, T_vb);

 [   [mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
    [mean(TP_cv),std(TP_cv),mean(TP_vb),std(TP_vb)]
    [mean(TN_cv),std(TN_cv),mean(TN_vb),std(TN_vb)]
    [mean(MSE_cv),std(MSE_cv),mean(MSE_vb),std(MSE_vb)]
]
end
result{kk}=[[mean(TP_cv),std(TP_cv);mean(TP_vb),std(TP_vb)], ...
    [mean(TN_cv),std(TN_cv);mean(TN_vb),std(TN_vb)],...
        [mean(MSE_cv),std(MSE_cv); mean(MSE_vb),std(MSE_vb)],...
    [mean(time_cv),std(time_cv);mean(time_vb),std(time_vb)]];
end
save(['SF_',Dyn,'_N',num2str(N),'_M',num2str(M),'.mat'],'result')
    
r=cell2mat(result');
%%
close all
clear
clc

N = 100; % Number of nodes
M = 200;
Dyn = 'C';
Alpha = -linspace(2,3,6);   % Alpha of the scale-free graph
maxiter = 100;
param.v0 = 1;
param.omega = 10^3;
param.r = 2+rand(N);

TP_cv = zeros(1,maxiter); time_cv = TP_cv; TN_cv = TP_cv; DegCorr_cv = TP_cv; 
TP_vb = zeros(1,maxiter); time_vb = TP_vb; TN_vb = TP_vb; DegCorr_vb = TP_vb; 
result = cell(1,10);
scale=0.1;
for kk=1:6
    exp_pow = Alpha(kk);
for numiter=1:maxiter
    display(uint16([N,M,scale,exp_pow,numiter]))

%define node degree distribution:
XAxis  = unique(round(logspace(0,log10(N),1000)));
YAxis  = unique(round(logspace(0,log10(N),1000))).^(exp_pow);
% create the graph with the required node degree distribution:
Graph = mexGraphCreateRandomGraph(N,XAxis,YAxis,0);
g=graph(Graph.Data(:,1),Graph.Data(:,2));
A=full(adjacency(g));
l=size(A,1);
if l<N
    A(1,(l+1):N)=1;
    A((l+1):N,1)=1;
end
    % simulate ele model 
    [F, O, W] = communication(A,M);
    F = F + randn(size(F))*scale;
    % Reconstruct network
    T_cv = zeros(size(A)); temp_cv = zeros(1,N);
    T_vb = zeros(size(A)); temp_vb = zeros(1,N); 
    for i=1:N
        % cv+lasso
        try
            tic;[B,info] = lasso(O,F(:,i),'CV',5);temp_cv(i)=toc;
        catch
            tic;[B,info] = lasso(O,F(:,i),'CV',5);temp_cv(i)=toc;
        end
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

    time_cv(numiter) = sum(temp_cv); time_vb(numiter) = sum(temp_vb); 
    [TP_cv(numiter), TN_cv(numiter), DegCorr_cv(numiter)] = Metrics(A, T_cv2);
    [TP_vb(numiter), TN_vb(numiter), DegCorr_vb(numiter)] = Metrics(A, T_vb);

 [   [mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
    [mean(TP_cv),std(TP_cv),mean(TP_vb),std(TP_vb)]
    [mean(TN_cv),std(TN_cv),mean(TN_vb),std(TN_vb)]
    [mean(MSE_cv),std(MSE_cv),mean(MSE_vb),std(MSE_vb)]
]
end
result{kk}=[[mean(TP_cv),std(TP_cv);mean(TP_vb),std(TP_vb)], ...
    [mean(TN_cv),std(TN_cv);mean(TN_vb),std(TN_vb)],...
        [mean(MSE_cv),std(MSE_cv); mean(MSE_vb),std(MSE_vb)],...
    [mean(time_cv),std(time_cv);mean(time_vb),std(time_vb)]];
end
save(['SF_',Dyn,'_N',num2str(N),'_M',num2str(M),'.mat'],'result')
    
r=cell2mat(result');
% m = MSE_cv;
% m(m>300) = mean(m(m<=300));
% [mean(m),std(m)]
% 
% close all
% clear
% clc
% 
% N = 50;
% M = 100;
% maxiter = 100;
% param.v0 = 1;
% param.omega = 10^3;
% param.r = 2+rand(N);
% scale = 0.1;
% exp_pow = 3;
% TP_cv = zeros(1,maxiter); time_cv = TP_cv; TN_cv = TP_cv; DegCorr_cv = TP_cv; 
% TP_vb = zeros(1,maxiter); time_vb = TP_vb; TN_vb = TP_vb; DegCorr_vb = TP_vb; 
% for numiter=1:maxiter
%     display([N,M,scale,exp_pow,numiter])
%     DegSeq = random('Generalized Pareto', 1./(exp_pow-1), 1, 0, [N,1]);
%     DegSeq = round(DegSeq); DegSeq(DegSeq<=0)=1; DegSeq(DegSeq>100)=100;
%     A = cm_net(DegSeq);
%     % simulate ele model 
%     [I, V] = ElectricalCurrent(A,M,param);
%     I = I + randn(size(I))*scale;
%     % Reconstruct network
%     T_cv = zeros(size(A)); temp_cv = zeros(1,N);
%     T_vb = zeros(size(A)); temp_vb = zeros(1,N); 
%     for i=1:N
%         G = V - V(:,i);
%         % cv+lasso
%         tic;[B,info] = lasso(G,I(:,i),'CV',10);temp_cv(i)=toc;
%         T_cv(i,:) = B(:,info.Index1SE);
%         % bayes recon
%         G(:,i)=[];
%         tic;model = vbr(G,I(:,i));temp_vb(i)=toc;
%         temp = model.theta';
%         T_vb(i,:) = [temp(1:(i-1)),0,temp(i:end)];
%     end
%     T_cv2=T_cv;T_cv2(T_cv2~=0)=1;
%     T_vb(T_vb>0.5)=1;T_vb(T_vb<=0.5)=0;
% 
%     time_cv(numiter) = sum(temp_cv); time_vb(numiter) = sum(temp_vb); 
%     [TP_cv(numiter), TN_cv(numiter), DegCorr_cv(numiter)] = Metrics(A, T_cv2);
%     [TP_vb(numiter), TN_vb(numiter), DegCorr_vb(numiter)] = Metrics(A, T_vb);
%     
%     [mean(time_cv),std(time_cv),mean(time_vb),std(time_vb)]
%     [mean(TP_cv),std(TP_cv),mean(TP_vb),std(TP_vb)]
%     [mean(TN_cv),std(TN_cv),mean(TN_vb),std(TN_vb)]
%     [mean(DegCorr_cv),std(DegCorr_cv),mean(DegCorr_vb),std(DegCorr_vb)]
% %     save(['N',num2str(N),'_M',num2str(M),'_m',num2str(m),'#',num2str(scale),'.mat'])
% end
% 
% [[mean(TP_cv),std(TP_cv);mean(TP_vb),std(TP_vb)], ...
%     [mean(TN_cv),std(TN_cv);mean(TN_vb),std(TN_vb)],...
%     [mean(time_cv),std(time_cv);mean(time_vb),std(time_vb)]]