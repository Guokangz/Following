hold on
clc
clear all

N = 2000;
T = 100;
Tf = 2*T;
t = linspace(0,Tf,N);
omega = (pi / T) ;
r = 1;
phi_0 = pi;
t_0 = 0;
p_x_val = p_x(r, omega, phi_0, t_0);
p_y_val = p_y(r, omega, phi_0, t_0);
alpha = cacu_alpha(p_x_val, p_y_val);
Phi_0 = atan(1 ./ (p_x_val + 1j * p_y_val));
psi0 =[-sin(Phi_0/2), cos(Phi_0/2)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,rho]=ode45('equation',t,psi0,opts,[],omega);
ita1_0 = [-sin(Phi_0/2), cos(Phi_0/2)];
ita2_0= [cos(Phi_0/2), sin(Phi_0/2)];
init_f = fidlity(ita1_0,psi0)
for nn = 1:N
t_0 = t(nn);
p_x_val = p_x(r, omega, phi_0, t_0);
p_y_val = p_y(r, omega, phi_0, t_0);
xdata(nn) =p_x_val;
ydata(nn) =p_y_val;
alpha = cacu_alpha(p_x_val, p_y_val);
Phi_0 = atan(1 ./ (p_x_val + 1j* p_y_val));
ita1 = [-sin(Phi_0/2), cos(Phi_0/2)];
ita2 = [cos(Phi_0/2), sin(Phi_0/2)];
rho_norm = rho(nn,:)/norm(rho(nn,:));
f1(nn) = fidlity(ita1,rho(nn,:));
f2(nn) = fidlity(ita2,rho(nn,:));
end
subplot(2,2,1)
plot(t/T,f1,t/T,f2,'LineWidth',2)

set(gca,'linewidth',1.5)
set(gca,'TickDir','in')
set(gca,'fontsize',20);
set(gca,'FontName','Times')
title('','fontsize',20,'interpreter','latex')
ylabel('$ |\langle \eta_k|\psi_1(t) \rangle|^2$','fontsize',20,'interpreter','latex')
xlabel('$\omega t / \pi$','fontsize',20,'interpreter','latex')
set(gca,'YLim',[0 1.01]);
set(gca,'XLim',[0 Tf/T]);
box on
grid on

function y = p_x(r,omega,phi_0,t)
    y = r*sin(omega*t+phi_0);
end
function y = p_y(r,omega,phi_0,t)
   y = 1-r*cos(omega*t+phi_0);
end
function c= cacu_alpha(x,y)
    
    Phi_c= atan(1./(x+1j*y));
    phi_r_c = real(Phi_c);
    phi_i_c = imag(Phi_c);
    alpha = x.*sinh(phi_i_c)./(sin(phi_r_c));
    c = alpha;
end 
function y = fidlity(state1,state2)
 state1_t = conj(state1);
 state2 = transpose(state2); 
 y = abs(dot(state1_t,state2))^2;
end 

