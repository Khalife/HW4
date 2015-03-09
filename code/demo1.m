
%% One iteration that updates \nabla f(\theta_n) using <GradSto>
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
   sigma = sqrt( rand ); % initial sigma_0
    beta = zeros(p,1);
    nnzeronumb = floor( .04*p );
    nnzeroentr = 8 * rand(nnzeronumb,1) + 1;
    nnzeroindx = unidrnd(p,nnzeronumb,1);% it is possible to get two equal indexes but it is very rare for large p: at least 1/p^2
    beta(nnzeroindx) = nnzeroentr;

theta = [beta; sigma]; % theta_0


lambda = 30;
gamma = 0.005;

%% GradSto : sample \nabla l(\theta)

for n = 1 : 10
Nm(n) = 200+n; 
Hnew = GradSto(Nm(n), theta, Z, X, Y);

% -- compute ERR, SEN and PRE (n, beta, BETA) ---
ERR(n) = norm( beta - BETA) / norm(BETA);

% -- proximal operator: P --
% min( -l(theta) + lambda |g|_1 : \nabla f(\theta_{n+1} = - Hnew .

 theta = P( theta + gamma*Hnew, gamma,lambda );
 beta = theta(1:end-1);


end

figure(); plot(ERR);
figure(); plot(BETA, '*-'); hold on;
plot( beta , 'ro-');



