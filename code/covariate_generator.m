function X=covariate_generator(N,p,rho)
%function covariate_generator(N,p,rho)
%   generates autoregressive covariate vectors
%   N is the vector dimension
%   p is the the number of vectors
%   rho parameter

    X=zeros(N,p);
    %% Initialization
    X(:,1)=randn(N,1);
    %% Auto-regression
    
    for i=2:p
        X(:,i)=rho*X(:,i-1)+sqrt(1-rho^2)*randn(N,1);
    end
    
end 