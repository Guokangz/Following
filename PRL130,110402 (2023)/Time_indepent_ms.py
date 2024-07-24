import time
import matplotlib.pyplot as plt
import numpy as np
from qutip import (Options, Qobj, basis, destroy, expect,mcsolve, mesolve, sigmax, sigmaz)

delta_max = 10 * 2 * np.pi

delta_min = -290* 2 * np.pi


t_list = np.linspace(0.0, 30, 1000)

T1 = 30

def delta_t(t,args=None):

    return  (delta_min+delta_max)/2-((delta_max-delta_min)/2)* np.cos((2*np.pi *t)/T1)


gamma2 = 880# dephasing rate
omega = 29  * 2 * np.pi
psi0 = basis(2, 0)  # initial state

sx = sigmax()
sz = sigmaz()

ee1 = [[0, 0], [0, 1]]

ee = Qobj(ee1)

ge1 = [[0, 1], [0, 0]]

ge = Qobj(ge1)




# perform the calculation for each combination of eps and A, store the result
# in a matrix


c_ops = np.sqrt(gamma2) * ge
# Define H1 (first time-dependent term)    # String method:

H0 = (omega/2.0)*sx

H1 = [ee, delta_t]

H = [H0, H1]

options = Options()
options.atol = 1e-6  # reduce accuracy to speed
options.rtol = 1e-5  # up the calculation a bit
options.nsteps = 100000000
options.rhs_reuse = True      
result0 = mesolve(H, psi0, t_list, c_ops, [ee],options = options)

print(c_ops)

plt.plot(t_list/T1, result0.expect[0], label="Expectation Value")
# 添加标题和轴标签
plt.title("Expectation Value vs Time")
plt.xlabel("Time")
plt.ylabel("Expectation Value")

# 添加图例
plt.legend()

# 显示图形
plt.show()
