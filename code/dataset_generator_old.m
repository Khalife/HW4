function [Y]=dataset_generator(N,p,q)
%function [Y]=dataset_generator(N,p,q)
%   generates a data set of random effect logistic regression observations
%   the regressor is sparce (98% of entries are set to zero)
%   the covariates are autoregressive with white noise
%   N is the dimension of an instance
%   p is the number of covariates
%   q is the dimension of loading vectors

    %% The covariate vectors
    rho=.8;
    X=covariate_generator(N,p,rho);
    
    %% The regressor
    beta=zeros(N,1);
    nnzeronumb=floor(.02*N);
    nnzeroentr=4*rand(nnzeronumb,1)+1;
    nnzeroindx=unidrnd(N,nnzeronumb,1);% it is possible to get two equal indexes but it is very rare for large p: at least 1/p^2
    beta(nnzeroindx)=nnzeroentr;
    
    %% The variance
    
    sigma=.1;
    
    %% The loading vector
    
    Z=zeros(q,p);
    
    for i=1:p
        indx=ceil(i*q/N);
        if indx<= q
            Z(indx,i)=1;
        end
    end
    
    %% The noise 
    
    U=randn(q,1);
    
    %% The dataset
    
    alphat= X'* beta+ sigma* Z' * U;
    
    logistic=@(x) exp(x)./(1+exp(x));
    
    alp=logistic(alphat);
    
    test=rand(p,1);
    
    Y=alp>test;
    
end