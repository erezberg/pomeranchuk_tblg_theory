% Compute surface of 1st order transition in nu, B, T phase space for simple thermodynamic model 
%% Parameters
e0 = 72;
dmu = 64;
dgamma = - (1/5)^2;



%% Nu, T
B = 0:1:20;


[nu, T] = meshgrid(1:0.02:1.5, 1e-2:0.02:(20+1e-2));

figure; hold on;
xlabel('\nu');
ylabel('T');

C = colormap;
for j = 1:length(B)
    DF = e0 - nu*dmu - T.*log(2*cosh(B(j)./T/2)) - 0.5*dgamma*T.^2;
    contour(nu, T, DF, [0,0], 'color', C(10*j,:),'linewidth', 1);
end

axis([1.0 1.15 0 14])
set(gca,'fontsize',16)


%% B, T 

[B,T] = meshgrid(0:0.2:14, 1e-2:0.02:(15+1e-2));
nu_dr = (e0 - T.*log(2*cosh(B./T/2)) - 0.5*dgamma*T.^2)/dmu;
figure; hold on; 
xlabel('B');
ylabel('T');
contour(B,T,nu_dr,20,'linewidth',1);
set(gca,'fontsize',16)

%% Mesh

[B,T] = meshgrid(0:1:12, 2:1:14);
nu_dr = (e0 - T.*log(2*cosh(B./T/2)) - 0.5*dgamma*T.^2)/dmu;
figure; hold on;
mesh(nu_dr,B,T,'linewidth',2);
xlabel('\nu');
ylabel('B');
zlabel('T', 'interpreter', 'latex');
set(gca,'fontsize',16)


