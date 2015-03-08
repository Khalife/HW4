function  Y = SampleG2(lambda, t )

% Y ~ g2 I_{y>=t}

% (1) sample from X ~ E(lambda)
X = exprnd(1/lambda);

% (2) t + X
% density of Y : lamda*exp(lamda*t) exp(-lamda*y) I_{y>=t} 
% it is a truncated exponential distribution
Y = t + X; 



