function Xg = SampleG( lambda, z, K, t, alpha)

mu = 1/z;

% F_IG = @(x,lambda,mu) cdf('normal',sqrt(lambda/x)*(x/mu -1), 0, 1) + ...
%         exp(2*lambda/mu)*cdf('normal',-sqrt(lambda/x)*(x/mu + 1) );
 
%     alpha = (exp(-K*t)/K) / ( exp(-K*t)/K + F_IG(t,lambda,mu)*4*exp(-z)/pi );
    
    U = rand;

    if U <= alpha
        Xg = SampleG2(K,t);
    else
        Xg = SampleG1(lambda,mu,t);
    end
