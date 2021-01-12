param.nf = 4;
param.W  = 1;   % Bandwidth (charge neutrality to top of band)
param.U  = 1.0; % Interaction strength
%param.E0 = 1;
param.N =   1000;   % number of discretized energies
param.eps_max = 4*param.W;
param.mu_max = 3*param.W;
param.nrand = 10; % Number of random initial guesses for each iteration
param.B = 0.04; % Zeeman field
param.Bv = param.B*[1 1 -1 -1];

param.eps = -param.eps_max:param.W/param.N:param.eps_max;

param.mu =  -param.mu_max:param.W/param.N:param.mu_max;

param.n = zeros(size(param.mu));
param.f = zeros(size(param.mu));

%% DOS
param.nu = abs(param.eps)/param.W^2;

%% Truncate DOS to bandwidth
param.nu = param.nu.*((param.eps<=param.W).*(param.eps>=-param.W));

%% Noralize nu such that n(W)=1
ntot = trapz(param.eps, param.nu);
param.nu = param.nu*2/ntot;


%% minimize free energy
phi =  0:0.02:5.5; % Chemical potentials
T = [4 6 8 10 12 14 16]/150*param.W;  % Temperatures. Below, the entropy is calculated and plotted at the even T's assuming they are equally spaced

mu_min = zeros(param.nf, length(phi),length(T));
n_min =  zeros(param.nf, length(phi),length(T));
F_min  = zeros(length(phi),length(T));

mu0 = [0.01,0.011,0.012,0.013]'; % Initial guess for mu

options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 3000, 'MaxIterations', 3000);

F_tr =   zeros(1,param.nrand+1);
mu_tr =  zeros(param.nf,param.nrand+1);

%% Main loop
tic
for jt = 1:length(T)
    fprintf('T = %.4f\n',T(jt))            
    param.T = T(jt);
    %% For each termperature, generate param.mu, param.n, and param.f
    [param.n,param.f] = init_param_finite_T(param);
    %% Loop over phi
    for jp = 1:length(phi)
        param.phi = phi(jp);
        %% Restrict the mu's to be in the range [phi-param.mu_max,phi+param.mu_max] 
        I = find(param.n==1,1);
        if (~isempty(I))    
            mu_max = param.mu(I);
        else
            mu_max = param.mu_max;
        end
        lb = - param.phi - mu_max*ones(4,1);
        ub = - param.phi + mu_max*ones(4,1);
        %% Minimize free energy
        F = @(x)F_mu(x,param);  % Define f = F_mu as anonymous function
        mu0_tr = - param.phi + 2*param.W*(rand(param.nf,param.nrand)-0.5);   % Random initial guess
        for jr = 1:(param.nrand+1)
            if (jr==1)
                mu_init = mu0;
            else
                mu_init = mu0_tr(:,jr-1);
            end
            [mu_tr(:,jr),F_tr(jr)]= fmincon(F,mu_init,[],[],[],[],lb,ub,[],options);
        end
        [~,I] = min(F_tr);
        mu_min(:,jp,jt) = mu_tr(:,I);
        F_min(jp,jt) = F_tr(I);
        mu0 = mu_min(:,jp,jt);
        n_min(:,jp,jt) = n_mu(mu_min(:,jp,jt),param);
        if (mod(jp+1,50)==0)
            fprintf('phi = %.4f\n',phi(jp))              
        end
    end
end
toc

n_min = sort(n_min,1);
n_tot = squeeze(sum(n_min,1));

%% Show entropy S = -dF/dT
figure; set(gca, 'fontsize', 14)
plot(n_tot(:,2:2:end), -diff(F_min(:,1:2:end),1,2)/(T(3)-T(1)),'linewidth',2);
xlabel('\nu');
ylabel('s');
set(gca, 'fontsize', 16)

%% Show chemical potential vs. density
figure; 
plot(n_tot, phi','linewidth',1)
xlabel('\nu')
ylabel('\phi')
set(gca, 'fontsize', 16)



