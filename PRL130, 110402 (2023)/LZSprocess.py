import matplotlib.pyplot as plt
import numpy as np
from qutip import (Options,Qobj,about, destroy, expect, num, propagator,
                   propagator_steadystate,basis, sigmax, sigmaz,mesolve)
# set up the parameters and start calculation
delta = 29 * 2 * np.pi  # qubit sigma_x coefficient


gamma2 = 800 # dephasing  rate

A_list= np.linspace(-290, 10,100)*2* np.pi

# pre-calculate the necessary operators
sx = sigmax()
sz = sigmaz()
sm = destroy(2)
sn = num(2)
custom_matrix_data1 = [[0, 0], [0, 1]]
custom_matrix_data2 = [[0, 1], [0, 0]]
ge=  Qobj(custom_matrix_data2)
ee = Qobj(custom_matrix_data1)
# initial state
psi0 = basis(2, 0)

tlist = np.linspace(0, 3000, 100)

# collapse operators: relaxation and dephasing
c_ops = [np.sqrt(gamma2) * ge]


# ODE settings (for list-str format)
options = Options()
options.atol = 1e-6  # reduce accuracy to speed
options.rtol = 1e-5  # up the calculation a bit
options.nsteps = 10000000

# perform the calculation for each combination of eps and A, store the result
# in a matrix

H0 = delta / 2.0 * sx

    # Define H1 (first time-dependent term)    # String method:
H1 = [ee/ 2, "A"]
    # Function method:
    # H1 = [- sz / 2, lambda t, args: args['eps'] ]
H = [H0, H1]
res_list = []

for A in A_list:
        
    args = {"A":A}

    res = mesolve(H, psi0, tlist, c_ops, [ee],args={"A": A},options = Options(nsteps = 10000000))
    
    res_list.append(res.expect[0]) 







fig, ax = plt.subplots(figsize=(8, 8))

plt.plot(tlist, res.expect[0], label="Expectation Value")

# 添加标题和轴标签
plt.title("Expectation Value vs Time")
plt.xlabel("Time")
plt.ylabel("Expectation Value")

# 添加图例
plt.legend()

# 显示图形
plt.show()


