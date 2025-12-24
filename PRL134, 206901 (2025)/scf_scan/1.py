import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
import io

# ==========================================
# 1. 物理常数与参考能量
# ==========================================
Ry_to_meV = 13605.698
# E_ref = E_slab + E_H2
E_ref = -2273.0948136022 

# ==========================================
# 2. 导入高密度扫描数据 (Flat Configuration)
# ==========================================
flat_data_str = """
2.000000   -2273.08070403
2.250000   -2273.08898638
2.500000   -2273.09447857
2.750000   -2273.09750157
3.000000   -2273.09876491
3.250000   -2273.09903449
3.500000   -2273.09883718
3.750000   -2273.09839376
4.000000   -2273.09783963
4.250000   -2273.09729052
4.500000   -2273.09680421
4.750000   -2273.09640242
5.000000   -2273.09607657
5.250000   -2273.09582024
5.500000   -2273.09562040
"""

# 读取数据
data = np.loadtxt(io.StringIO(flat_data_str))
dist_flat = data[:, 0]  # 距离 (A)
E_tot_flat = data[:, 1] # 总能 (Ry)

# 计算吸附能 (meV)
E_ads_flat = (E_tot_flat - E_ref) * Ry_to_meV

# ==========================================
# 3. 绘图设置
# ==========================================
fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

# --- 曲线平滑处理 ---
X_smooth = np.linspace(dist_flat.min(), dist_flat.max(), 300)
spl = make_interp_spline(dist_flat, E_ads_flat, k=3)
Y_smooth = spl(X_smooth)

# --- 绘制主曲线 ---
# 1. 红色实线 (拟合线)
ax.plot(X_smooth, Y_smooth, 'r-', linewidth=3, label='H$_2$ on Ag(111) (Flat)')
# 2. 原始数据点 (红色空心圆)
ax.plot(dist_flat, E_ads_flat, 'o', color='red', markersize=8, 
        markerfacecolor='white', markeredgewidth=2, zorder=5)

# --- 关键数据标注 ---
# 找到最低点
min_idx = np.argmin(E_ads_flat)
min_dist = dist_flat[min_idx]
min_E = E_ads_flat[min_idx]

# 标注最低点数值
ax.annotate(f'Min: {min_E:.1f} meV\n@ {min_dist:.2f} $\AA$', 
            xy=(min_dist, min_E), xytext=(min_dist+0.5, min_E+20),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=2),
            fontsize=16, color='red', fontweight='bold')

# ==========================================
# 4. 样式美化
# ==========================================
# 0能级参考线
ax.axhline(0, color='gray', linestyle=':', linewidth=1.5) 

# 标签
ax.set_xlabel(r'Distance from Surface ($\AA$)', fontsize=18, fontweight='bold')
ax.set_ylabel(r'Adsorption Energy (meV)', fontsize=18, fontweight='bold')

# 设置范围 (稍微收紧一点，突出势阱)
ax.set_xlim(1.9, 5.6)
ax.set_ylim(-80, 150) # 上限设高一点展示排斥区

# 刻度设置 (朝内，大号字体)
ax.tick_params(direction='in', which='major', length=8, width=1.5, labelsize=14, top=True, right=True)
ax.tick_params(direction='in', which='minor', length=4, width=1.0, top=True, right=True)

# 图例
ax.legend(frameon=False, fontsize=14, loc='upper right')

# 网格
ax.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()

# 保存图片
# plt.savefig('PES_Flat_Clean.png', dpi=300)
plt.show()