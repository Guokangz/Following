function rhodot=masterequation(t,rho,~, dim,J,L,N)

x=reshape(rho,dim,1);

%% 先不考虑相位的事情
%% H = Omega_1*srr*sD' + Omega_2*sD*s11' + (Omega_1*srr*sD' + Omega_2*sD*s11')' + Delta_2*s11*s11'; %% 原始

for i = 1:N 

    if mod(i, 2) == 0
    H(i,i+1) = J*(0.1+0.8*(1-exp(-3*(L-t)/L)/(1-exp(-3))));
    H(i+1,i) = J*(0.1+0.8*(1-exp(-3*(L-t)/L)/(1-exp(-3))));
    else
    H(i,i+1) = J*(0.1+0.8*(1-exp(-3*t/L)/(1-exp(-3))));
    H(i+1,i) = J*(0.1+0.8*(1-exp(-3*t/L)/(1-exp(-3))));
    end

end

a1=-1i*H*x;

rho1=reshape(a1,dim,1);
rhodot=rho1;


