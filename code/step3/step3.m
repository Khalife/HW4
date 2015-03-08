% Step 3 code
addpath('hw1');
s=@(x) exp(x)/(1+exp(x));

N=100;
q=4;
p=5;
z=zeros(N,q);
sigma=0.5;
beta=zeros(p,1);
X=zeros(N,p);
Y=zeros(N,1);

%init 
K_chaine=1000; %longueur chaine 
w=ones(N,K);


for k=2:K_chaine
    G=zeros(q);
    m=zeros(q,1);
    for j=1:N
        G=G+w(j,k-1)*transpose(z(j,:))*z(j,:);
        m=m+((Y(j)-1/2)-w(j,k-1)*dot(X(j,:),beta)*transpose(z(j)));
    end    
    Gamma=inv(eye(q)+sigma^2*G);
    mu=sigma*Gamma*m;
    u=mvnrnd(mu,Gamma);
    
    for i=1:N
        w(i,k)=(1/4)*HomeWork1(abs(dot(X(i,:),beta)+sigma*dot(z(i,:),u)));
    end
    k
end

%GradSto
v=zeros(p+1,N);
gradL=zeros(p+1,1);
for i=1:N
    v(:,i)=[transpose(X(i,:));dot(z(i,:),u)];
    gradL=gradL+v(:,i)*(Y(i)-s(dot(X(i,:),beta)+sigma*dot(z(i,:),u)));
end

gradL=(1/N)*gradL; %gradient en theta = [beta, sigma];
