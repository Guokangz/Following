#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import  pi, sqrt
from pulse import Pulse
import tqdm
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
mpl.rcParams['font.size'] = 20
N = 100
tf_list =np.linspace(115,1500,N)
gamma = 0.3*pi

#%%
# 计算脉冲
pulse1 = []
pulse2 = []
pulse3 = []
for i in tqdm.tqdm(range(N)):

    tf = tf_list[i]*1e-9
    delta1 = 1*pi*1e9
    pulse = Pulse(tf,gamma,delta1,None)
    pulse1.append(max(pulse.Omega1())/1e7*pi)
    delta1 = 1.5*pi*1e9
    pulse = Pulse(tf,gamma,delta1,None)
    pulse2.append(max(pulse.Omega1())/1e7*pi)
    delta1 = 2*pi*1e9
    pulse = Pulse(tf,gamma,delta1,None)
    pulse3.append(max(pulse.Omega1())/1e7*pi)
#%%
# 绘制脉冲图

plt.figure(figsize=(6, 12))
plt.plot(tf_list,pulse1,linestyle='dotted',color= 'green',label = r'$\mathrm{\Delta = 1\pi GHz}$')
plt.plot(tf_list,pulse2,linestyle='--',color= 'b',label = r'$\mathrm{\Delta = 1.5\pi GHz}$')
plt.plot(tf_list,pulse3,linestyle='-',color= 'r',label = r'$\mathrm{\Delta = 2\pi GHz}$')
plt.xlabel(r'$\mathrm{t_{f}(ns)}$')
plt.ylabel(r'Amplitude of $\widetilde{\Omega}_{1,4}(\mathrm{\pi MHz)}$'),
plt.grid(True)
plt.legend(framealpha=0.1) 
plt.tight_layout() 
plt.savefig('fig5a.png', dpi=600, format='png')
plt.show()


# %%
