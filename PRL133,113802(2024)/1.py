import numpy as np
import matplotlib.pyplot as plt
def cacu_alpha(x,y):
    Phi = np.arctan(1/(x+1j*y))
    phi_r = np.real(Phi)
    phi_i = np.imag(Phi)
    alpha = x*np.sinh(phi_i)/np.sin(phi_r)
    return alpha
# 定义 x 和 y 的范围
x_range = np.linspace(-1.5, 1.5, num=200)
y_range = np.linspace(-2.5, 2.5, num=200)

# 初始化存储 alpha 值的二维数组
alpha_values = np.zeros((len(x_range), len(y_range)))

# 双循环遍历 x 和 y 的范围
for i, x_val in enumerate(x_range):
    for j, y_val in enumerate(y_range):
        alpha_values[j, i] = cacu_alpha(x_val, y_val)

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt

# 创建一个新的图形对象
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 创建网格
X, Y = np.meshgrid(x_range, y_range)
Z = alpha_values

# 绘制三维图
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.plot_surface(X, Y, -Z, cmap='viridis')
# 设置标签
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Alpha Values')

# 显示图形
plt.show()