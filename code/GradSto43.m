
function  gradL = GradSto43(w,u, theta, Z, X, Y, opt)

sigma = theta(end); 
beta = theta(1:end-1);

s = opt.s; 

auxv = @( X,Y, Z, u ) Y - s(X*beta + sigma*Z*u);
% input Nx1, Nxp, Nxq , qx1 ---> output : Nx1

[N , p] = size( X );
q = size(Z, 2);


%init
K_chaine = size( w, 2);
% w = ones(N, K_chaine);
% u = zeros( q, K_chaine);

Hk = zeros( p + 1 ,K_chaine );

for k = 1 : K_chaine
    % calculate H(uk) = sum_i auxv_i (x_i ; z'*u)
    tmp = [X'; (Z*u(:,k))']; % first p components// w.r.t \beta
    auxk = auxv( X,Y,Z, u(:,k) );
    tmp = repmat( auxk', p+1 ,1).*tmp;
    
    Hk(:,k) =  sum( tmp, 2 );
    
end


gradL = mean(Hk(:,end-100:end) , 2 );



