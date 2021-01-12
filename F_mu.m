function F = F_mu(mu,param)
%% Compute F = F_var(\mu+\phi) + U/2 n(\mu+\phi) n(\mu+\phi) + \mu n(\mu+\phi) 
F = 0;
%% Interpolate (param.mu, param.f) for each flavour
if (abs(param.B)<1e-10)
    for j=1:param.nf
        if (abs(mu(j)+param.phi)<param.mu_max)
            F = F + interp1(param.mu,param.f,mu(j)+param.phi);
        else
            % Note that here f is assumed to be a symmetric function
            F = F + param.f(1) - (abs(mu(j)+param.phi)-param.mu_max);
        end
    end
    %% Interpolate n
    n = zeros(param.nf,1);
    for j=1:param.nf
        if (abs(mu(j)+param.phi)<param.mu_max)
            n(j) = interp1(param.mu,param.n,mu(j)+param.phi);
        else
            n(j) = sign(mu(j)+param.phi);
        end    
    end
else
    for j=1:param.nf
        phi_j = param.phi+param.Bv(j)/2;
        if (abs(mu(j)+phi_j)<param.mu_max)
            F = F + interp1(param.mu,param.f,mu(j)+phi_j);
        else
            % Note that here f is assumed to be a symmetric function
            F = F + param.f(1) - (abs(mu(j)+phi_j)-param.mu_max);
        end
    end
    %% Interpolate n
    n = zeros(param.nf,1);
    for j=1:param.nf
        if (abs(mu(j)+phi_j)<param.mu_max)
            n(j) = interp1(param.mu,param.n,mu(j)+phi_j);
        else
            n(j) = sign(mu(j)+phi_j);
        end    
    end    
end
    
%% Add interaction and mu parts to F
F = F + 0.5*param.U*(n'*(ones(param.nf,1)*sum(n)-n));
F = F + mu'*n;
end
