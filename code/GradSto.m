
function Hnew = GradSto(Nmax, theta, Z, X, Y, opt)
% Step 3 code

sigma = opt.SIGMA; % theta(end); % for 4.3
beta = theta;   % (1:end-1); % for 4.3

s = @(x) exp(x)./(1+exp(x));

auxv = @( X,Y, Z, u ) Y - s(X*beta + sigma*Z*u); 
% input Nx1, Nxp, Nxq , qx1 ---> output : Nx1 

[N , p] = size( X );
q = size(Z, 2);


%init 
K_chaine = Nmax; %1000; %longueur chaine 
w = ones(N, K_chaine);

Hk = zeros( p ,1 );

for k = 2:K_chaine
    G=zeros(q);
    m=zeros(q,1);
    for j=1:N
        G = G+w(j,k-1)*transpose(Z(j,:))*Z(j,:);
        m = m + ((Y(j)-1/2) - w(j,k-1)* dot(X(j,:),beta)) *transpose(Z(j,:));
    end    
    Gamma = inv(eye(q)+sigma^2*G);
    mu = sigma*(eye(q)+sigma^2*G)\m;
    u = mvnrnd(mu,Gamma)';
    
    for i = 1:N
        w(i,k) = 4*HomeWork1( 0.5*abs(dot(X(i,:),beta)+sigma*dot(Z(i,:),u )));
    end
    disp(k);
% calculate H(uk) = sum_i auxv_i (x_i ; z'*u)
tmp = X'; 
auxk = auxv( X,Y,Z, u );
tmp = repmat( auxk', p,1).*tmp;

Hk = Hk + sum( tmp, 2 );

end

%GradSto

Hnew = Hk / (K_chaine - 1 );


