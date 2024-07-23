import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import symbols, Matrix,simplify

# 定义符号变量
r_eff, omiga, delta = sp.symbols('r_eff omiga delta')

# 定义刘维尔量矩阵
H1 = sp.Matrix([[-r_eff, 1j*omiga/2, -1j*omiga/2, 0],
                [1j*omiga/2, -((r_eff/2)+1j*delta), 0, -1j*omiga/2],
                [-1j*omiga/2, 0, -((r_eff/2)-1j*delta), 1j*omiga/2],
                [r_eff, -1j*omiga/2, 1j*omiga/2, 0]])


# 求解本征值
eigenvalues_H1 = H1.eigenvals()
simplified_eigenvalues_H1 = [eigenvalue.simplify() for eigenvalue in eigenvalues_H1]

delta_values = np.linspace(-300, 300, num=50, endpoint=True)
gamma_values = np.linspace(50, 850, num=50, endpoint=True)

E1_data = np.zeros((len(delta_values), len(gamma_values)))
E2_data = np.zeros((len(delta_values), len(gamma_values)))


for i, dita in enumerate(delta_values):

    for j, gamma in enumerate(gamma_values):
        # 创建变量值字典

        variable_sd = {r_eff: gamma, delta:dita, omiga:29 * 2 * np.pi}
        variables = {r_eff: gamma, delta:dita * 2 * np.pi, omiga:29 * 2 * np.pi}
        print(variable_sd)
        # 替换符号值并计算本征值
        E1 = simplified_eigenvalues_H1[1].subs(variables)
        E2 = simplified_eigenvalues_H1[2].subs(variables)
        
        # 使用 evalf() 方法将符号表达式转换为数值\
        E1_numeric = E1.evalf()
        E2_numeric = E2.evalf()
        print(E2_numeric)
        print(E1_numeric)
        E1_i =sp.im(E1_numeric)
        E2_i =sp.im(E2_numeric)
        print(E1_i)
        print(E2_i)
        E1_data[i, j] = E1_i
        E2_data[i, j] = E2_i
    


delta_grid, gamma_grid = np.meshgrid(gamma_values,delta_values)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 绘制曲面图
surface1 = ax.plot_surface(gamma_grid, delta_grid, E1_data, cmap='viridis',rcount=200, ccount=200)
surface2 = ax.plot_surface(gamma_grid, delta_grid, E2_data, cmap='plasma',rcount=200, ccount=200)

# 添加颜色条
fig.colorbar(surface1, ax=ax, label='E1')
fig.colorbar(surface2, ax=ax, label='E2')


ax.set_xlabel('$\u0394$/2$\u03c0$(kHz)')
ax.set_ylabel('$\u03b3$(kHz)')
ax.set_zlabel('Energy(kHz)')

plt.show()