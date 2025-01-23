
clc
clear all
% r = 0.5;
% phi_0 = pi;
 M = 2;
N =100;
% omega = (pi / T) ;
% x = r*sin(omega*t+phi_0);
% y = 1-r*cos(omega*t+phi_0);
% phi_c= atan(1/(x+1j*y));
% phi_r = real(phi_c);
% phi_i = imag(phi_c);
% alpha = x*sinh(phi_i)/(sin(phi_r))
% ;

x_list = linspace(-1.5,1.5,N);
y_list = linspace(-3,3,N);
H = zeros(M,M);

for xx = 1:N
    x_tem = x_list(xx);
    for yy = 1:N
        y_tem = y_list(yy);
        phi_c= atan(1/(x_tem+1i*y_tem));
        phi_r = real(phi_c);
        phi_i = imag(phi_c);
        alpha = x_tem*sinh(phi_i)/(sin(phi_r));
        H(1,2) = alpha*sin(phi_c);
        H(2,1) = alpha*sin(phi_c);
        H(1,1) = alpha*cos(phi_c);
        H(2,2) = alpha*-cos(phi_c);
        [V,D] = eig(H);
        data1(yy,xx) = real(D(1,1));
        data2(yy,xx) = real(D(2,2));
    end
end
surface(x_list,y_list,data1)
hold on 
surface(x_list,y_list,data2)
shading interp
grid on
