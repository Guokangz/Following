#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, pi, sqrt
from tqdm import tqdm
from pulse import Pulse
# 生成需要的hamiltonian space
v1 = np.array([0, 0, 0, 0, 0])
v2 = np.array([0, 0, 0, 0, 0])
# H_delta
delta = 2*pi*1e9
H_delta = Qobj(np.diag([0,delta,0,delta,0]))
# H_Omega1
matrix = np.outer(v1, v2)
matrix[[1,0],[0,1]] = 1
H_Omega1 = Qobj(matrix)
# H_Omega2
matrix = np.outer(v1, v2)
matrix[[2,1],[1,2]] = 1
H_Omega2 = Qobj(matrix)
# H_Omega3
matrix = np.outer(v1, v2)
matrix[[3,2],[2,3]] = 1
H_Omega3 = Qobj(matrix)
# H_Omega4
matrix = np.outer(v1, v2)
matrix[[4,3],[3,4]] = 1
H_Omega4 = Qobj(matrix)


cop_matrix = np.ones((5, 5))
gamma1 = 0.01e3
gamma2 = 30e3
gamma3 = 0.06e3
gamma4 = 30e3
gamma5 = 0
gamma_a = np.array([gamma1,gamma2,gamma3,gamma4,gamma5])

for i in range(5):

    for j in range(5):

        cop_matrix[i, j] = (gamma_a[i] +gamma_a[j])/2

c_ops = Qobj(sqrt(cop_matrix))

psi0 = basis(5, 0)
psif = basis(5, 4)

N =100
Nstep = 1000
tf = 1000e-9
t = tf
gamma_1 = np.linspace(0.001, 0.5, N)
times = np.linspace(0, t, 10000)
all_probabilities2 = []
fid = []
for i in tqdm(range(N)):
    opts = Options(nsteps=100000, atol=1e-6)
    gamma = gamma_1[i]*pi
    pulse = Pulse(tf,gamma,delta, None)
    args = { 'tf': tf,'delta':delta,'gamma':gamma}
    H2 = [H_delta,[H_Omega1,pulse.Omega1()],
     [H_Omega2,pulse.Omega2()],
     [H_Omega3,pulse.Omega3()],
     [H_Omega4,pulse.Omega4()]] 
    # with the counterdiabatic driving
    result2 = mesolve(H2, psi0, times, args=args,options=opts)
    fid.append(fidelity(psif, result2.states[-1]))

#%%
plt.plot(gamma_1,fid)
plt.xlabel(r'$\Omega0$')
plt.ylabel('Probability')
plt.grid(True)
plt.show()


 # %%
