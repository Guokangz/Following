#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt
import matplotlib as mpl
from tqdm import tqdm
from pulse import Pulse
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
mpl.rcParams['font.size'] = 20
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
# H_Omegacd
matrix = np.outer(v1, v2)
matrix[4,0] = -1
matrix[0,4] = 1
H_Omegacd = Qobj(matrix)

cop_matrix = np.ones((5, 5))
gamma1 = 0.01e3
gamma2 = 30e3
gamma3 = 0.06e3
gamma4 = 30e3
gamma5 = 0e3
gamma_a = np.array([gamma1,gamma2,gamma3,gamma4,gamma5])

for i in range(5):

    for j in range(5):

        cop_matrix[i, j] = (gamma_a[i] +gamma_a[j])/2

c_ops = Qobj(sqrt(cop_matrix))

H_Omegacd = 1j * H_Omegacd
psi0 = basis(5, 0)
N =50
tlist = np.linspace(100, 1000, N)
all_probabilities1 = []
all_probabilities2 = []
for i in tqdm(range(N)):
   
    tf = tlist[i]*1e-9
    sigma = tf/6
    tao = tf/10
    times = np.linspace(0, tf, 1000)    
    omega0 = 30*pi*1e6
    opts = Options(nsteps=100000, atol=1e-6)
    pulse = Pulse(tf, omega0, tao, sigma, None)
    args = {'tao': tao, 'tf': tf, 'sigma': sigma, 'omega0': omega0}
    H1 = [H_delta,[H_Omega1,pulse.Omega1()],
     [H_Omega2,pulse.Omega2()],[H_Omega3,pulse.Omega3()],
     [H_Omega4,pulse.Omega4()],[H_Omegacd,pulse.Omegacd()]]
    result1 = mesolve(H1, psi0, times, c_ops, args=args,options=opts)
    omega0 = 125*pi*1e6
    pulse = Pulse(tf, omega0, tao, sigma, None)
    H2 = [H_delta,[H_Omega1,pulse.Omega1()],
     [H_Omega2,pulse.Omega2()],
     [H_Omega3,pulse.Omega3()],
     [H_Omega4,pulse.Omega4()]] 
    # with the counterdiabatic driving
    result2 = mesolve(H2, psi0, times, c_ops, args=args,options=opts)
    states = [basis(5, i) for i in range(5)]

    projectors = [psi * psi.dag() for psi in states]
    probabilities1 = np.array([expect(P, result1.states) for P in projectors])
    probabilities2 = np.array([expect(P, result2.states) for P in projectors])
    all_probabilities1.append(probabilities1[-1][-1])
    all_probabilities2.append(probabilities2[-1][-1])

#%%

plt.figure(figsize=(12, 6))
plt.plot(tlist, all_probabilities1,linewidth=2, color='red',label='C-STIRSAP')
plt.plot(tlist, all_probabilities2,'--',linewidth=2, color='blue',label='C-STIRAP')
plt.xlabel(r'$\mathrm{t_f}$(ns)')
plt.ylabel('Transfer efficiency')
plt.grid(True)
plt.legend() 
plt.tight_layout() 
plt.savefig('fig3b.png', dpi=600, format='png')
plt.show()


# %%
