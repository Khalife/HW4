function X = SampleG1(lambda, mu, t)

 % inverse of psi_1 for algo3
Ipsi = @(y,mu) mu + 0.5*mu^2 *y - 0.5*mu*sqrt(4*mu*y + (mu*y).^2);

if mu > t
    % call algorithm 2
    while true
        while true
            E = exprnd(1/lambda);
            Eb = exprnd(1/lambda);
            if E^2 <= 2*Eb/t;
                break;
            end
        end
        X = t/(1+t*E)^2;
        alpha = exp(-0.5*X/mu^2);
        U = rand;
        if U <= alpha
            break;
        end
    end
    
else
    % call algorithm 3
    while true
        Y = randn;
        Y = Y^2;
        X = Ipsi(Y,mu);
        U = rand;
        if U > mu/(mu+X)
            X = mu^2/X;
        end
        if X <= t
            break;
        end
    end
end
    
    
