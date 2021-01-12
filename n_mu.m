% Calculate n as a function of mu(:) by interpolation
function n = n_mu(mu,param)
    n = zeros(param.nf,1);
    for j=1:param.nf
        if (abs(mu(j)+param.phi)<param.mu_max)
            n(j) = interp1(param.mu,param.n,mu(j)+param.phi);
        else
            n(j) = sign(mu(j)+param.phi);
        end    
    end
end
