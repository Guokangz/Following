import numpy as np
from scipy.integrate import ode

# 定义微分方程组函数
def diff_eqs(t, rho, omega, gamma_eff, delta):

    rho11, rho12, rho21, rho22 = rho

    drho11_dt = 1/2 * (1j * omega * (rho12 - rho21) - np.abs(gamma_eff) * (rho11 - 2 * rho22))
    drho12_dt = 1j * delta * rho12 + 1/2 * 1j * omega * (rho11 - rho22)
    drho22_dt = 1/2 * (1j * omega * (-rho12 + rho21) - np.abs(gamma_eff) * rho22)
    drho21_dt = -1j * delta * rho21 - np.abs(gamma_eff) * rho21 + 1/2 * 1j * omega * (-rho11 + rho22)

    return [drho11_dt, drho12_dt, drho21_dt, drho22_dt]

omega = 29*2*np.pi
gamma_eff = 800
rho0 = [1, 0, 0, 0]

# 定义时间点
t_start = 0
t_end = 30
num_points = 100
t_values = np.linspace(t_start, t_end, num_points)

def delta_value(t):

    return (-290 + 10 * t) * 2 * np.pi
# 创建ode对象


# 初始化 results 列表
results = []

# 数值求解微分方程组
for t in t_values:

    delta = delta_value(t)  # 获取 delta(t) 的值
    print(delta)
    r = ode(diff_eqs).set_integrator('zvode', method='bdf')
    r.set_initial_value(rho0, t_start).set_f_params(omega, gamma_eff, delta)
    result = r.integrate(t)
    results.append(result)
    
    
    while r.successful() and r.t < t:
        result = r.integrate(t)
        results.append(result)

# 解析结果
results = np.array(results)
rho22_result = results[:, 3]

# 打印结果或绘制图形
print(rho22_result)



