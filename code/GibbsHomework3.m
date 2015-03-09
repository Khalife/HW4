
function  [w, u] = GibbsHomework3(Nmax, theta, Z, X, Y, opt)

% stop = 0;
if opt.sigmaknown
sigma = opt.SIGMA;
beta = theta;
else
    sigma = theta(end);
    beta = theta(1:end-1);
end
    

% s = @(x) exp(x)./(1+exp(x));
% auxv = @( X,Y, Z, u ) Y - s(X*beta + sigma*Z*u);
% % input Nx1, Nxp, Nxq , qx1 ---> output : Nx1

[N , p] = size( X );
q = size(Z, 2);


%init
K_chaine = Nmax; %1000; %longueur chaine
w = ones(N, K_chaine);
u = zeros( q, K_chaine);

for k = 2:K_chaine
% find G, m with equivalent compact form
    G=zeros(q);
    m=zeros(q,1);
%             for j=1:N
%                 G = G+w(j,k-1)*transpose(Z(j,:))*Z(j,:);
%                 m = m+ ( (Y(j)-1/2)-w(j,k-1)* dot(X(j,:),beta) ) *transpose(Z(j,:)) ;
%             end

    wZ = repmat(w(:,k-1),1,q).*Z;
    G = Z'*wZ;
    m = sum( repmat( ((Y - 0.5) - w(:,k-1).*(X*beta)), 1,q).* Z, 1)';
% end of finding G, m
    
    Gamma = inv(eye(q)+sigma^2*G);
    mu = sigma * (eye(q)+sigma^2*G)\ m;
    u(:,k-1) = mvnrnd(mu,Gamma)';
    
    for i = 1:N
        w(i,k) = 4*HomeWork1( 0.5*abs( dot(X(i,:),beta)+sigma*dot(Z(i,:),u(:,k-1) )));
    end
    disp(k);
        
end

    wZ = repmat(w(:,end),1,q).*Z;
    G = Z'*wZ;
    m = sum( repmat( ((Y - 0.5) - w(:, end).*(X*beta)), 1,q).* Z, 1)';
    
    Gamma = inv(eye(q)+sigma^2*G);
    mu = sigma * (eye(q)+sigma^2*G)\ m;
    u(:,end) = mvnrnd(mu,Gamma)';

% after this : GradSto

% Hnew = mean(Hk(:,end-100:end) , 2) ;% / (K_chaine - 1 );


