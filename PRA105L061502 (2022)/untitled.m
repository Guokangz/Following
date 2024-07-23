eta = 11.76;
gamma = -0.335;

d = linspace(6,16,100);
J = eta * exp(gamma*d);

plot(d,J)