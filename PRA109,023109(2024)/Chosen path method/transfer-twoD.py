#%%
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, pi, sqrt
import numdifftools as nd
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

""" cop_matrix = np.ones((5, 5))

gamma1 = 0.01
gamma2 = 30
gamma3 = 0.06
gamma4 = 30
gamma5 = 0
gamma_a = np.array([gamma1,gamma2,gamma3,gamma4,gamma5])

for i in range(5):

    for j in range(5):

        cop_matrix[i, j] = (gamma_a[i] +gamma_a[j])/2

c_ops = Qobj(cop_matrix) """

psi0 = basis(5, 0)
N =25
tlist = np.linspace(100, 1500, N)
gamma_list = np.linspace(0.001, 0.5, N)
all_probabilities = np.zeros((N, N))
#%%
for i in tqdm(range(N)):
   
   for j in tqdm(range(N)):
        opts = Options(nsteps=10000, atol=1e-6)
        gamma = gamma_list[j]*pi
        tf =  tlist[i]*1e-9
        times = np.linspace(0,tf,1000)
        pulse = Pulse(tf,gamma,delta, None)
        args = { 'tf': tf, 'gamma': gamma}
        H2 = [H_delta,[H_Omega1,pulse.Omega1()],
        [H_Omega2,pulse.Omega2()],
        [H_Omega3,pulse.Omega3()],
        [H_Omega4,pulse.Omega4()]] 
        # with the counterdiabatic driving
        result2 = mesolve(H2, psi0, times,args=args,options=opts)
        states = [basis(5, i) for i in range(5)]
        projectors = [psi * psi.dag() for psi in states]
        #probabilities1 = np.array([expect(P, result1.states) for P in projectors])
        probabilities2 = np.array([expect(P, result2.states) for P in projectors])
        #all_probabilities1.append(probabilities1[-1][-1])
        all_probabilities[-j,i] = probabilities2[-1][-1]

#%%
plt.pcolormesh(all_probabilities, cmap='rainbow',vmin=0.95, vmax=1)
#plt.contourf(all_probabilities, cmap='rainbow', levels=5)
# 添加颜色条
plt.colorbar()

# 显示图像
plt.show()
# %%
plt.imshow(all_probabilities, cmap='rainbow', vmin=0.96, vmax=1)

# 绘制等高线图
contours = plt.contour(all_probabilities, colors='black')
plt.clabel(contours, inline=True, fontsize=8)

# 添加颜色条
plt.colorbar()

# 显示图像
plt.show()
# %%
