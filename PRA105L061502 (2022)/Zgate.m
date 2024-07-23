clc; clear all; close all;
tic
% ÉèÖÃËã×Ó

dim = 11;

%% ÑÝ»¯Ëã·û
%% ³õÌ¬ ºÍ Ä©Ì¬

zero = zeros(dim,1);
zero(1,1) =1;
psi_0= zero;
zero = zeros(dim,1);
zero(11,1) =1;
psi_f = zero;

%%
L = 1;
N = 10;
J = 1;
time = 10;
T = linspace(0,time,1000);
p0 = psi_0;
options = odeset('reltol',1e-6,'abstol',1e-6);

[t,rho] = ode45('masterequation',T,p0,options,dim,J,L,N);

for n = 1:1000

    pideal = reshape(rho(n,:),dim,1);
    Fidelity(n) = abs(psi_f'*pideal)^2;
    n;
end

legend
hold on
plot(T/time,Fidelity,'LineWidth',2);
legend
grid minor
toc


