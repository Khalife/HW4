function w = HomeWork1(z)

% Shuyu DONG - HW1 Simulation based learning
% shuyu.dong@polytechnique.edu

%% Parameters and functions
lambda = 1;
mu = 1/z;
K = pi^2/8 + z^2/2;

t = 0.64;

a = @(n,x,t) pi*(n+1/2)*(2/(pi*x))^1.5 * exp( -(2*(n+1/2)^2)/x )*(x<=t)*(x>0)...
           + pi*(n+1/2)* exp(- ((n+1/2)^2*pi^2*x) / 2 )*(x>t);

Sn0 = @(x,z,t) cosh(z)*exp(-x*z^2/2) * a(0,x,t);

% CDF of the Inverse Gaussian
F_IG = @(x,lambda,mu) cdf('normal',sqrt(lambda/x)*(x/mu -1), 0, 1) + ...
        exp(2*lambda/mu)*cdf('normal',-sqrt(lambda/x)*(x/mu + 1) );    
% rejection parameter
alpha = (exp(-K*t)/K) / ( exp(-K*t)/K + F_IG(t,lambda,mu)*4*exp(-z)/pi );


%% Algo 4 : generation of w ~ f_*

while true
n = 0;
% sample from g
X = SampleG(lambda, z, K, t, alpha);

S0 = Sn0(X,z,t);
U = S0 *rand;
while true
    n = n + 1;
    Sn = S0 + (-1)^n*a(n,X,t);
    if ( U<=Sn && ~mod(n,2) ) || (U>Sn && mod(n,2))
        break;
    end
end

if ~ mod(n,2)
break;
end

end

w = X;



