
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
%% 
opt.sigma = sqrt( rand ); % initial sigma_0
    beta = zeros(p,1);
    nnzeronumb = floor( .04*p );
    nnzeroentr = 8 * rand(nnzeronumb,1) + 1;
    nnzeroindx = unidrnd(p,nnzeronumb,1);% it is possible to get two equal indexes but it is very rare for large p: at least 1/p^2
    beta(nnzeroindx) = nnzeroentr;

%% GradSto : sample \nabla l(\theta)

theta = [beta; opt.sigma]; % theta_0

for n = 1 : 1
Nm(n) = 200+n; 
Hnew = GradSto(Nm(n), theta, Z, X, Y);
end
