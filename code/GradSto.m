function Hnew = GradSto(Nmax, theta, z, X, Y)
% Step 3 code

sigma = theta(end);
beta = theta(1:end-1);

s = @(x) exp(x)./(1+exp(x));

auxv = @( X,Y, Z, u ) Y - s(X*beta + sigma*Z*u); 
% input Nx1, Nxp, Nxq , qx1 ---> output : Nx1 

[N , p] = size( X );
q = size(z, 2);


%init 
K_chaine = Nmax; %1000; %longueur chaine 
w = ones(N, K_chaine);
u = zeros(q, K_chaine);

Hk = zeros( p+1 ,1 );

for k = 2:K_chaine
    G=zeros(q);
    m=zeros(q,1);
    for j=1:N
        G=G+w(j,k-1)*transpose(z(j,:))*z(j,:);
        m=m+((Y(j)-1/2)-w(j,k-1)* dot(X(j,:),beta) *transpose(z(j,:)));
    end    
    Gamma = inv(eye(q)+sigma^2*G);
    mu = sigma*Gamma*m;
    u(:,k-1) = mvnrnd(mu,Gamma);
    
    for i = 1:N
        w(i,k)=(1/4)*HomeWork1( 0.5*abs(dot(X(i,:),beta)+sigma*dot(z(i,:),u(:,k-1) )));
    end
    disp(k);
% calculate H(uk) = sum_i auxv_i (x_i ; z'*u)
tmp = [X'; (z*u(:,k-1))']; 
auxk = auxv( X,Y,z,u(:,k-1) );
tmp = repmat( auxk', p+1,1).*tmp;

Hk = Hk + sum( tmp, 2 );

end

%GradSto

Hnew = Hk / (K_chaine - 1 );


