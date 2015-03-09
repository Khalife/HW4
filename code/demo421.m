
%%  4.2 - Strategie 2 - Algo (5) - \nabla f(\theta_n) using <GradSto>
%
close all;
clearvars except N;
if ~ exist('N','var')
addpath( './step3/hw1');
end
N = 500;
p = 1000;
q = 5;

%% dataset and ture parameters
SIGMA = sqrt( 0.1 );
opt.SIGMA = SIGMA;
[Y, X, Z, BETA] = dataset_generator(N,p,q, opt);
Bnz = (BETA ~=0 );
Bsumnz = sum( Bnz);
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

% theta_0: IN 4.2 sigma is assumed to be known
% theta = [beta; sigma];
theta = beta;
opt.sigmaknown = 1;
%% hyperparameters
% iterations for algo (5)
n = 20;

gamma = 0.005;
Nm = 150 + (1:n) ; % may lift to 200+ later...
lambdas = 30; %[10 30 90 ];

sen = @( b ) sum( (b~=0).* Bnz )/Bsumnz ;
pre = @( b ) sum( (b~=0).* Bnz ) /sum( (b~=0) );
ERR = zeros(1,n);
SEN = zeros(1,n);
PRE = zeros(1,n);
s = @(x) exp(x)./(1+exp(x));
opt.s = s;

for i = 1 : length( lambdas )
lambda = lambdas(i);

%% GradSto : sample \nabla l(\theta)

for t = 1 : n
    disp( ['algo5 iteration: t = ',int2str(t)] );
%     gradL  = GradSto(Nm(t), theta, Z, X, Y, opt);
    [w, u ] = GibbsHomework3(Nm(t), theta, Z, X, Y, opt);
    gradL = GradSto(w,u, theta, Z, X, Y, opt);
%     if stop
%         break;
%     end
    % -- compute ERR, SEN and PRE (n, beta, BETA) ---
    ERR(t) = norm( beta - BETA) / norm(BETA);
    SEN(t) = sen( beta );
    PRE(t) = pre( beta );
    % -- proximal operator: P --
    % min( -l(theta) + lambda |g|_1 : \nabla f(\theta_{n+1} = - gradL .
    
    theta = P( theta + gamma*gradL, gamma,lambda );
    beta = theta ; % In 4.2, sigma is known
end

figure(); plot(ERR);
figure(); plot(SEN);
figure(); plot(PRE);

figure(); plot(BETA, '*-'); hold on;
plot( beta , 'ro-');
pause(.01);
end

savestep4 = 'step4/figures/';
saveas(1, [savestep4,'421-ERRnew.jpg'] );
saveas(2, [savestep4,'421-SENnew.jpg'] );
saveas(3, [savestep4,'421-PREnew.jpg'] );
saveas(4, [savestep4,'421-betaBETAnew.jpg'] );


