
%%  4.2 - Strategie 2 - Algo (5) - \nabla f(\theta_n) using <GradSto>
%
addpath( './step3/hw1');
N = 500;
p = 1000;
q = 5;

%% dataset and ture parameters
SIGMA = sqrt( 0.1 );
opt.SIGMA = SIGMA;
[Y, X, Z, BETA] = dataset_generator(N,p,q, opt);

%% Proximal P

P = @ (u, gamma, lambda) (u-gamma*lambda).*(u>= gamma*lambda) + ...
    (u+gamma*lambda).*(u<= -gamma*lambda);

%% Initial parameters:
% initial sigma
sigma = sqrt( rand );
% initial beta:
beta = zeros(p,1);
nnzeronumb = floor( .04*p );
nnzeroentr = 8 * rand(nnzeronumb,1) + 1;
nnzeroindx = unidrnd(p,nnzeronumb,1);% it is possible to get two equal indexes but it is very rare for large p: at least 1/p^2
beta(nnzeroindx) = nnzeroentr;

% theta_0:
theta = [beta; sigma];

%% hyperparameters
% iterations for algo (5)
n = 10;

gamma = 0.01 ./ sqrt(1:n);
Nm = 100 + ceil( sqrt(1:n) ); % may lift to 200+ later...
lambdas = 30; %[10 30 90 ];

for i = 1 : length( lambdas )
lambda = lambdas(i);

%% GradSto : sample \nabla l(\theta)

for t = 1 : n
    
    Hnew = GradSto(Nm(t), theta, Z, X, Y);
    
    % -- compute ERR, SEN and PRE (n, beta, BETA) ---
    ERR(t) = norm( beta - BETA) / norm(BETA);
    
    % -- proximal operator: P --
    % min( -l(theta) + lambda |g|_1 : \nabla f(\theta_{n+1} = - Hnew .
    
    theta = P( theta + gamma(t)*Hnew, gamma(t),lambda );
    beta = theta(1:end-1);
end

figure(); plot(ERR);
figure(); plot(BETA, '*-'); hold on;
plot( beta , 'ro-');
pause(.01);
end


