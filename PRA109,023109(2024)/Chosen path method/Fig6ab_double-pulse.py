
#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import  pi, sqrt
from pulse import Pulse
from pulse2 import Pulse2

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

Nstep = 10000
tf = 1000e-9
t = tf
times  = np.linspace(0,t,Nstep)
gamma = 0.3*pi
#%%
import matplotlib.pyplot as plt
# 计算脉冲
pulse1 = Pulse(tf,gamma,delta)
pulse2 = Pulse2(tf,gamma,delta)

pulse1_1 = pulse1.Omega1()
pulse1_2 = pulse1.Omega2()
pulse1_3 = pulse1.Omega3()
pulse1_4 = pulse1.Omega4()

pulse2_1 = pulse2.Omega1()
pulse2_2 = pulse2.Omega2()
pulse2_3 = pulse2.Omega3()
pulse2_4 = pulse2.Omega4()
empty_array = np.zeros(100)
# 绘制脉冲图
pulse_1 = np.concatenate((pulse1_1, empty_array, pulse2_1))
pulse_2 = np.concatenate((pulse1_2, empty_array, pulse2_2))
pulse_3 = np.concatenate((pulse1_3, empty_array, pulse2_3))
pulse_4 = np.concatenate((pulse1_4, empty_array, pulse2_4))

x = np.linspace(0, 2*t, 20100)
plt.figure(figsize=(12, 6))
plt.plot(x*1e9,pulse_1/pi/1e6,color= 'green',linestyle = 'solid',label = r'$\mathrm{\widetilde{\Omega}_{1,4}}$')
plt.plot(x*1e9,pulse_2/pi/1e6,color= 'blue',linestyle = '-.',label = r'$\mathrm{\widetilde{\Omega}{2}}$')
plt.plot(x*1e9,pulse_3/pi/1e6,color= 'red',linestyle = '-.',label = r'$\mathrm{\widetilde{\Omega}_{3}}$')
plt.xlabel(r'$\mathrm{t\quad(ns)}$')
plt.ylabel(r'$\mathrm{Rabi\quadFrequency(\pi\quadMHz)}$')
plt.grid(True) 
plt.legend(framealpha=0.1)

plt.tight_layout() 
plt.savefig('fig6a.png', dpi=600, format='png')
plt.show()
#%%
print(x)
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
times = np.linspace(0, 2*t, 20100)
args = { 'tf': tf, 'delta': delta,'gamma':gamma}
#——布局图——————————————————————————————————————————————————————————————
opts = Options(nsteps=100000, atol=1e-6)
# without the counterdiabatic driving
H1 = [H_delta,[H_Omega1,pulse_1],
     [H_Omega2,pulse_2],
     [H_Omega3,pulse_3],
     [H_Omega4,pulse_4]]
result1 = mesolve(H1, psi0, times,c_ops, args=args,options=opts)
result2 = mesolve(H1, psi0, times, args=args,options=opts)
states = [basis(5, i) for i in range(5)]
projectors = [psi * psi.dag() for psi in states]
# 计算每个状态的概率
#%%
probabilities1 = np.array([expect(P, result1.states) for P in projectors])
probabilities2 = np.array([expect(P, result2.states) for P in projectors])

# 绘制概率随时间的变化
plt.figure(figsize=(12, 6))
colors = ['b', 'g', 'g', 'k', 'y']
for i, prob in enumerate(probabilities1):
    if i == 0:
        plt.plot(times*1e9, probabilities2[0],label=fr'$|{i+1}\rangle$''without decay',color = 'red',linestyle='--')
        plt.plot(times*1e9, prob, label=fr'$|{i+1}\rangle$''with decay',color = 'blue',linestyle='solid')
    else:
        plt.plot(times*1e9, prob, label=fr'$|{i+1}\rangle$', color=colors[i % len(colors)],linestyle='-.')
plt.xlabel('t(ns)')
plt.ylabel('Probability')
plt.legend(framealpha= 0.1,prop={'size': 13})
plt.grid(True)
plt.tight_layout() 
plt.savefig('fig6b.png', dpi=600, format='png')
plt.show()


# %%
