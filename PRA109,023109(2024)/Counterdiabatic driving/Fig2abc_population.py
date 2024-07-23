#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import  pi, sqrt
from pulse import Pulse
import matplotlib as mpl

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
# 创建哈密顿量


# %%

""" #构建脉冲函数
def Omega1(t, args):
    tao = args['tao']
    tf = args['tf']
    sigma = args['sigma']
    return sqrt(Omega2(t, args)**2 + Omega3(t, args)**2)  

def Omega2(t, args):
    tao = args['tao']
    tf = args['tf']
    sigma = args['sigma']
    return omega0 * np.exp(-((t-tao-tf/2)/sigma)**2)

def Omega3(t, args):
    tao = args['tao']
    tf = args['tf']
    sigma = args['sigma']
    return omega0 * np.exp(-((t+tao-tf/2)/sigma)**2)

def Omega4(t, args):
    tao = args['tao']
    tf = args['tf']
    sigma = args['sigma']
    return sqrt(Omega2(t, args)**2 + Omega3(t, args)**2) 

def Omega2_t(t):
    return Omega2(t, args)

def Omega3_t(t):
    return Omega3(t, args)

# 定义导数函数
def dOmega2_dt(t):
    return nd.Derivative(Omega2_t,step=1e-9)(t)

def dOmega3_dt(t):
    return nd.Derivative(Omega3_t,step=1e-9)(t)

def Omegacd(t, args):
    tao = args['tao']
    tf = args['tf']
    sigma = args['sigma']
    
    return (dOmega2_dt(t)*Omega3(t, args) - dOmega3_dt(t)*Omega2(t, args))/(Omega2(t, args)**2 + Omega3(t, args)**2)
 """

omega0 = 30*pi*1e6

tf = 1000e-9

sigma = tf/6

tao = tf/10

t = np.linspace(0, tf, 1000)

args = {'tao': tao, 'tf': tf, 'sigma': sigma}
# 计算脉冲
Pulse = Pulse(tf, omega0, tao, sigma, None)
pulse1 = Pulse.Omega1()/omega0
pulse2 = Pulse.Omega2()/omega0
pulse3 = Pulse.Omega3()/omega0
pulse4 = Pulse.Omega4()/omega0 
pulsecd = Pulse.Omegacd()/omega0 
# 绘制脉冲图
#%%
plt.figure(figsize=(12, 6))
plt.plot(t*1e9, pulse1, linestyle='-.', color='k')
plt.plot(t*1e9, pulse2, linestyle='-.', color='r', label=r'$\Omega_{2}/\Omega_{0}$')
plt.plot(t*1e9, pulse3, linestyle='-.', color='b', label=r'$\Omega_{3}/\Omega_{0}$')
plt.plot(t*1e9, pulse4, linestyle='-.', color='k', label=r'$\Omega_{1,4}/\Omega_{0}$')
plt.plot(t*1e9, pulsecd, linestyle='-', color='green', label=r'$\Omega_{\mathrm{cd}}/\Omega_{0}$')
plt.xlabel('t(ns)')
plt.ylabel('Rabi frequency')
plt.grid(True)  
plt.legend(framealpha=0.1)
plt.tight_layout() 
plt.savefig('fig2a.png', dpi=600, format='png')
plt.show()

# %%
# 合成哈密顿量
H_Omegacd = 1j * H_Omegacd
H = [H_delta,[H_Omega1,Pulse.Omega1()],
     [H_Omega2,Pulse.Omega2()],[H_Omega3,Pulse.Omega3()],
     [H_Omega4,Pulse.Omega4()],[H_Omegacd,Pulse.Omegacd()]]

# %%
# 定义初始状态
psi0 = basis(5, 0)

# 定义时间列表
times = np.linspace(0, tf, 1000)

# 定义耗散算子
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

#%%
times = np.linspace(0, tf, 1000)
#——布局图——————————————————————————————————————————————————————————————
opts = Options(nsteps=10000, atol=1e-6)
H1 = [H_delta,[H_Omega1,Pulse.Omega1()],
     [H_Omega2,Pulse.Omega2()],[H_Omega3,Pulse.Omega3()],
     [H_Omega4,Pulse.Omega4()],[H_Omegacd,Pulse.Omegacd()]]
result1 = mesolve(H1, psi0, times, args=args,options=opts)
# without the counterdiabatic driving
H2 = [H_delta,[H_Omega1,Pulse.Omega1()],
     [H_Omega2,Pulse.Omega2()],
     [H_Omega3,Pulse.Omega3()],
     [H_Omega4,Pulse.Omega4()]]
result2 = mesolve(H2, psi0, times, args=args,options=opts)
states = [basis(5, i) for i in range(5)]
projectors = [psi * psi.dag() for psi in states]
# 计算每个状态的概率
#%%
probabilities = np.array([expect(P, result1.states) for P in projectors])
colors = ['b', 'g', 'b', 'k', 'y']  # 创建颜色列表
for i, prob in enumerate(probabilities):
    if i == 4:
        plt.plot(times*1e9, prob,linestyle='-', color='r',label=fr'$|{i+1}\rangle$',)
    else:
        plt.plot(times*1e9, prob,linestyle='-.', color=colors[i % len(colors)],label=fr'$|{i+1}\rangle$')
plt.xlabel(r'$\mathrm{t{}(ns)}$')
plt.ylabel('Population')
plt.legend(framealpha=0.1)
plt.grid(True)
plt.tight_layout() 
plt.savefig('fig2b.png', dpi=600, format='png')
plt.show()



# %%
probabilities = np.array([expect(P, result2.states) for P in projectors])
colors = ['b', 'g', 'b', 'k', 'y']  # 创建颜色列表
for i, prob in enumerate(probabilities):
    if i == 4:
        plt.plot(times*1e9, prob*1e3,linestyle='-', color='r',label=fr'$|{i+1}\rangle\times10^3$',)
    else:
        plt.plot(times*1e9, prob,linestyle='-.', color=colors[i % len(colors)],label=fr'$|{i+1}\rangle$')
plt.xlabel(r'$\mathrm{t{}(ns)}$')
plt.ylabel('Population')
plt.legend(framealpha=0.1)
plt.grid(True)
plt.tight_layout() 
plt.savefig('fig2c.png', dpi=600, format='png')
plt.show()
# %%
