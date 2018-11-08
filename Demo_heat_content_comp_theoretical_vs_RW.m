%
%This demo program is used to show
% (1) how to generate one power law graph using Barabasi-Albert model;
% (2) compuate the theoretical heat content;
% (3) give a random walk simulation result of the heat content;
% (4) compare theoretical and random walk results.

% Output figure: heat contents (black: theoretical; red: random walk)

%Initialization
% N: Total number of nodes in the graph;
% mean_deg: mean degree of the graph;
% B_percentage: the percentage of boundary nodes;
% Time: heat content display time length;
% Interval: heat content display interval;
% Initial_walker: intial random walkers on each node;

N = 1000;
mean_deg = 10;
B_percentage = 0.02;
Time = 5; 
Interval = 0.01; 
Initial_walker = 10;

%-------------------------------------------------------------------------
%-------------------------- generate a B-A graph--------------------------

% Generate BA graph
% Output A_ba, the Adjacency matrix;
A_ba = BAModel(N, mean_deg);

% Randomly selecting the nodes with the smallest degrees as the boundaries. 
% Output: B_ba
B_ba = Boundary(B_percentage, A_ba);

% Compute eigenvalues and corresponding alpha value;
% Output: lambda_ba, alpha_ba;
[ lambda_ba, alpha_ba, D_ba] = coef_computing(A_ba, B_ba);



%-------------------------------------------------------------------------
%-------------------------- Heat content comparison ----------------------
disp('Displaying the heat content comparisons between theoretical and random walk...')
figure;
hold on
% Compute the theoretical heat content.
% Output: hc_ba
[ hc_ba ] = hc_Theoretical(lambda_ba,alpha_ba,Time,Interval,'k');
% Compute the random walk heat content.
% Output: hc_sim
[ hc_sim ] = hc_RWSimulation(A_ba, B_ba, Time, Interval, Initial_walker,'r');
hold off
legend('Theoretical','Random Walk')