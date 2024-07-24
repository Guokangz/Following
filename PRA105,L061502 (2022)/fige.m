T=100;
step = 1/T;

for m = 1:T

t = m*step; 
t
L = 1;
N = 10;
J = 1;


for i = 1:N 

    if mod(i, 2) == 0
    H(i,i+1) = J2(t,L,J);
    H(i+1,i) = J2(t,L,J);
    else
    H(i,i+1) = J1(t,L,J);
    H(i+1,i) = J1(t,L,J);
    end

end

[V, D] = eig(H);

mat(:,m)=diag(D);

end



function [result] =J1(t,L,J)

    result = J*(0.1+0.8*(1-exp(-3*t/L)/(1-exp(-3))));
    
end
function [result] =J2(t,L,J)

    result = J*(0.1+0.8*(1-exp(-3*(L-t)/L)/(1-exp(-3))));
    
end