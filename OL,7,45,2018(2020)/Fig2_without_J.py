#%%
import guan
import numpy as np
import matplotlib.pyplot as plt
from numpy import  pi, sqrt,cos
import matplotlib as mpl
import tqdm
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
mpl.rcParams['font.size'] = 20
# 生成需要的hamiltonian space
t =1
def H_ssh(N,omega):

    hamiltonian = np.zeros((N, N))

    for i in range(N):

        hamiltonian[i, i] = 0
       
    for i in range(N-1):
        hamiltonian[i, i+1] = (1+(-1)**i*cos(omega*t))
        hamiltonian[i+1, i] = (1+(-1)**i*cos(omega*t))

    import guan
    guan.statistics_of_guan_package()
    return hamiltonian

N = 5
start = 0
end =2
step = 0.01  # 设置步长，根据需要"调整
values = np.arange(start, end + step, step)
eigenvalues = []
for phistep in tqdm.tqdm(values):

    omega = phistep*pi
    hamiltonian = H_ssh(N,omega)
    eigvals  = guan.calculate_eigenvalue(hamiltonian)
    eigenvalues.append(eigvals) 

eigenvalues = np.array(eigenvalues).T
# %%
plt.figure(figsize=(8, 6)) 
for i in range(len(eigenvalues)):
    if i == 2:
        plt.plot(values, eigenvalues[i],'m',label='No.20')
    else:
        plt.plot(values, eigenvalues[i],'#808080')

plt.xlabel(r'$\theta/\pi$')
plt.ylabel('$\mathrm{E}$')
plt.grid(True)
plt.tight_layout() 
plt.savefig('fig2a.eps', dpi=600, format='eps')
plt.show()


# %%
P1=[]
S1=[]
start = 0
end =2
step = 0.001  # 设置步长，根据需要调整
values = np.arange(start, end + step, step)
for value in tqdm.tqdm(values):
    phi = value*pi
    hamiltonian = H_ssh(N,phi);
    eigenvectors =[]
    eigenvectors = guan.calculate_eigenvector(hamiltonian) 
    eigenvectors =eigenvectors.T
    eigenvectors =eigenvectors[2]
    P1.append(abs(eigenvectors))
P1 = np.array(P1).T
#%%
P1 = np.flip(P1)
plt.figure(figsize=(8, 6)) 
plt.pcolormesh(P1,cmap='hot')
plt.xlabel(r'$\theta/\pi$')
plt.ylabel(r'$\mathrm{Lattice\quadSite}$')
x_values = np.linspace(0, P1.shape[1],5)  # 创建 x 轴的刻度值
x_labels = [f'{np.round(value/1e3,1)}' for value in x_values]  # 创建 x 轴的刻度标签
plt.xticks(x_values, x_labels)  # 设置 x 轴的刻度和标签

y_values =[0.5,1.5,2.5,3.5,4.5]  # 创建 y 轴的刻度值
y_labels = [1,2,3,4,5]  # 创建 y 轴的刻度标签
plt.yticks(y_values, y_labels)
plt.colorbar()
plt.tight_layout() 
plt.savefig('fig2b.eps', dpi=600, format='eps')
plt.show()
# %%
