function [Y, X, Z, beta] = dataset_generator(N,p,q ,opt)
%function [Y]=dataset_generator(N,p,q)
%   generates a data set of random effect logistic regression observations
%   the regressor is sparce (98% of entries are set to zero)
%   the covariates are autoregressive with white noise
%   N is the dimension of an instance
%   p is the number of covariates
%   q is the dimension of loading vectors

    %% The covariate vectors
    rho = .8;
    X = covariate_generator(N,p,rho);
    
    %% The regressor
    beta = zeros( p,1 ); %N,1);
    nnzeronumb = floor(.02*p ); %N);
    nnzeroentr = 4*rand(nnzeronumb,1)+1;
    nnzeroindx = unidrnd(p ,nnzeronumb,1);% it is possible to get two equal indexes but it is very rare for large p: at least 1/p^2
    beta(nnzeroindx) = nnzeroentr;
    
    %% The variance
    if nargin < 4
    sigma = sqrt( .1 );
    else
        sigma = opt.SIGMA;
    end
    %% The loading vector
    
    Z = zeros(N, q); %zeros(q,p);
    
    for i = 1 : N %p
        indx = ceil(i*q/N);
         %if indx<= q
            Z(i, indx) = 1; %indx,i)=1;
         %end
    end
    
    %% The noise 
    
    U=randn(q,1);
    
    %% The dataset
    
    alphat= X* beta+ sigma* Z * U; % X'*beta +...
    
    logistic=@(x) exp(x)./(1+exp(x));
    
    alp=logistic(alphat);
    
    test=rand(N,1); %p,1);
    
    Y=alp>test;
    
end