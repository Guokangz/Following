#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tqdm as tqdm
from qutip import parallel_map
from qutip.ui.progressbar import TextProgressBar
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
mpl.rcParams['font.size'] = 20

def H_ssh():
    N = 5
    hamiltonian = np.zeros((N, N))

    hamiltonian[0, 1] = 1
    hamiltonian[1, 0] = 1
    hamiltonian[2, 3] = 1
    hamiltonian[3, 2] = 1
    return hamiltonian

def H_ssh2():
    N = 5
    hamiltonian = np.zeros((N, N))
    hamiltonian[1, 2] = 1
    hamiltonian[2, 1] = 1
    hamiltonian[3, 4] = 1
    hamiltonian[4, 3] = 1
    return hamiltonian

def H_J():
    N = 5
    hamiltonian = np.zeros((N, N))
    hamiltonian[0, 2] = 1
    hamiltonian[2, 0] = 1
    hamiltonian[2, 4] = 1
    hamiltonian[4, 2] = 1
    return hamiltonian



def H1_coeff(t,args):
    omega = args['omega']
    return 1 - np.cos(omega * t)

def H2_coeff(t,args):
    omega = args['omega']
    return 1 + np.cos(omega * t)

def HJ_coeff(J, args):
    J = args['J']
    return J


H_1 = Qobj(H_ssh())
H_2 = Qobj(H_ssh2())
HJ = Qobj(H_J())
H = [[H_1, H1_coeff], [H_2, H2_coeff], [HJ, HJ_coeff]]
print(H)
#%%
N = 200
fid=np.zeros((N,N))
J_list = np.linspace(0.1, 6, N)
omega_list = np.linspace(0.001, 0.2, N)

# %%
# %%

def compute_fidelity(args):
    H = args['H']
    J = args['J']
    omega = args['omega']
    t = np.pi / omega
    psi0 = basis(5,0)
    psif = basis(5,4)
    times = np.linspace(0, t, 1000)
    options = Options(nsteps=10000)
    result = sesolve(H, psi0, times, args={'t':t, 'omega':omega, 'J':J}, options=options)
    return fidelity(psif, result.states[-1])

args_list = [{'H': H, 'J': J, 'omega': omega} for J in J_list for omega in omega_list]
progress_bar = TextProgressBar(len(args_list))
results = parallel_map(compute_fidelity, args_list,progress_bar=progress_bar)
fid = np.array(results).reshape((N, N))
# %%
plt.figure(figsize=(10, 8))
plt.pcolormesh(fid.T, cmap='hot')
plt.xlabel(r'$J$')
plt.ylabel(r'$\Omega$')
plt.colorbar()
plt.title('Fidelity')
xticks = np.linspace(0, N, 5)  # 5 个刻度
xlabels = np.round(np.linspace(J_list[0], J_list[-1], 5))  # 对应的标签
plt.xticks(xticks, xlabels)
# 设置 y 轴的刻度值和标签
yticks = np.linspace(0, N, 5)  # 5 个刻度
ylabels = np.round(np.linspace(omega_list[0], omega_list[-1], 5),2)  # 对应的标签
plt.yticks(yticks, ylabels)
plt.tight_layout() 
plt.savefig('fig5a.eps', dpi=600, format='eps')
plt.show()
# %%
