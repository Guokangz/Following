#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import  pi, sqrt
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
# 创建哈密顿量

# %%
Nstep = 10000
tf = 1000e-9
t = tf
times  = np.linspace(0,t,Nstep)
gamma = 0.6*pi
#%%
import matplotlib.pyplot as plt
# 计算脉冲
pulse = Pulse(tf,gamma,delta)
pulse1 = pulse.Omega1()
pulse2 = pulse.Omega2()
pulse3 = pulse.Omega3()
pulse4 = pulse.Omega4()
# 绘制脉冲图
plt.figure()
plt.plot(times*1e9,pulse1/pi/1e6,color= 'green',linestyle = 'solid',label = r'$\mathrm{\widetilde{\Omega}_{1,4}}$')
plt.plot(times*1e9,pulse2/pi/1e6,color= 'red',linestyle = '-.',label = r'$\mathrm{\widetilde{\Omega}_{2}}$')
plt.plot(times*1e9,pulse3/pi/1e6,color= 'blue',label = r'$\mathrm{\widetilde{\Omega}_{3}}$',linestyle = '-.',)
plt.xlabel(r'$\mathrm{t (ns)}$')
plt.ylabel(r'$\mathrm{Rabi\quadFrequency (\pi MHz)}$')
plt.grid(True)  
plt.legend(framealpha=0.1)
plt.tight_layout() 
plt.savefig('fig5b.png', dpi=600, format='png')
plt.show()

# %%
# 定义初始状态
psi0 = basis(5, 0)

# 定义耗散算子
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
#c_ops = []

print(c_ops)
#%%
times = np.linspace(0, t, 10000)
args = { 'tf': tf, 'delta': delta,'gamma':gamma}
#——布局图——————————————————————————————————————————————————————————————
opts = Options(nsteps=100000, atol=1e-6)
# without the counterdiabatic driving
H1 = [H_delta,[H_Omega1,pulse1],
     [H_Omega2,pulse2],
     [H_Omega3,pulse3],
     [H_Omega4,pulse4]]
result1 = mesolve(H1, psi0, times,c_ops, args=args,options=opts)
states = [basis(5, i) for i in range(5)]
projectors = [psi * psi.dag() for psi in states]
# 计算每个状态的概率
#%%
probabilities = np.array([expect(P, result1.states) for P in projectors])
colors = ['b', 'g', 'g', 'k', 'y']
# 绘制概率随时间的变化
plt.figure()
for i, prob in enumerate(probabilities):
    if i == 4:
        plt.plot(times*1e9, prob, label=fr'$|{i+1}\rangle$', color = 'red',linestyle='--')
    else:
        plt.plot(times*1e9, prob,color=colors[i % len(colors)], label=fr'$|{i+1}\rangle$')
plt.xlabel(r'$\mathrm{t(ns)}$')
plt.ylabel('Population')
plt.legend(framealpha=0.1)
plt.grid(True)
plt.tight_layout() 
plt.savefig('fig5c.png', dpi=600, format='png')
plt.show()

# %%
