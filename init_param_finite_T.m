% Initialize parameters: param.n, param.f 
function [n,f] = init_param_finite_T(param)
    T = param.T;
    thetaf = param.eps<0;  % theta function Theta(-eps)
    thetaf((length(param.eps)-1)/2+1)=0.5; % Theta(0)=0.5

    %%Given grid of param.mu, generate param.n 
    n = zeros(1,length(param.mu));
    for j = 1:length(param.mu)
        n(j) = trapz(param.eps, param.nu.* ...
            (1./(1 + exp((param.eps - param.mu(j))/T)) - thetaf));
    end
    
    %Given grid of param.mu, generate param.f
    f = zeros(1,length(param.mu));    
    for j = 1:length(param.mu)
        exp_arg  = -(param.eps - param.mu(j))/T; % argument of exponential
        log_func = (exp_arg<=100).*log(1 + exp(min(exp_arg,100)))...
            + (exp_arg>100).*exp_arg;
        f(j) = -T*trapz(param.eps, param.nu.*(log_func ...
                              + (param.eps - param.mu(j))/T.*thetaf) );
    end    
end
