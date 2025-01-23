function rhodot = equation(t,rho,opts,flag,omega)
r = 1;
phi_0 = pi;
M = 2;
x = r*sin(omega*t+phi_0);
y = 1-r*cos(omega*t+phi_0);
phi_c= atan(1/(x+1j*y));
phi_r = real(phi_c);
phi_i = imag(phi_c);
alpha = x*sinh(phi_i)/(sin(phi_r));

H = zeros(M,M);
H(1,2) = sin(phi_c);
H(2,1) = sin(phi_c);
H(1,1) = cos(phi_c);
H(2,2) = -cos(phi_c);

H = alpha*H;

 a1=-1i*H*rho;

rhodot=a1;
