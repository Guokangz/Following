#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import  pi
import tqdm as tqdm

rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
mpl.rcParams['font.size'] = 20
v1 = np.array([0, 0, 0, 0, 0])
v2 = np.array([0, 0, 0, 0, 0])
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
def HJ_coeff(t,args):
    J = args['J']
    return J


H_1 = Qobj(H_ssh())
H_2 = Qobj(H_ssh2())
HJ = Qobj(H_J())
H = [[H_1, H1_coeff], [H_2, H2_coeff], [HJ, HJ_coeff]]
print(H)
#%%
J = 6
omega = 5e-6
t = pi/omega
psi0 = basis(5,1)
psif = basis(5,3)
times = np.linspace(0, t, 10000)
args ={'t':t, 'omega':omega, 'J':J}
options = Options(nsteps=10000)
result = sesolve(H, psi0, times,args=args, options=options)

v1 = np.array([0, 0, 0, 0, 0])
for i in result.states:
    v1 = np.vstack((v1, i.full().T))
v2 = np.flip(abs(v1))
v2_transposed = np.transpose(v2)
#%%
plt.figure(figsize=(10, 8))
plt.pcolormesh(v2_transposed, cmap='hot')
plt.xlabel(r'$\Omega t$')
plt.ylabel(r'$\mathrm{Lattice\quadSite}$')
x_values = np.linspace(0,v2_transposed.shape[1],3)  # 创建 x 轴的刻度值
x_labels = [0,0.5,1]  # 创建 x 轴的刻度标签
plt.xticks(x_values, x_labels)  # 设置 x 轴的刻度和标签

y_values =[0.5,1.5,2.5,3.5,4.5]  # 创建 y 轴的刻度值
y_labels = [1,2,3,4,5]  # 创建 y 轴的刻度标签
plt.yticks(y_values, y_labels)
plt.colorbar()
plt.tight_layout() 
#plt.savefig('fig5b.eps', dpi=600, format='eps')
plt.show()
# %%

