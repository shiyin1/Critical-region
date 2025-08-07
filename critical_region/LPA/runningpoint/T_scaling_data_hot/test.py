import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# 整合所有数据点
all_x = [0.1, 0.2, 0.3, 0.4]  # 横坐标共有4组
all_y = [               # 每组的纵坐标长度不同
    np.array([1, 2, 3]),
    np.array([1.5, 2.5, 3.5, 4.5]),
    np.array([2, 4, 6]),
    np.array([1, 3, 5, 7, 9])
]
all_z = [               # 每组的测量值，长度跟对应的 y 相同
    np.array([10, 20, 15]),
    np.array([12, 22, 18, 16]),
    np.array([14, 24, 19]),
    np.array([13, 23, 17, 15, 11])
]

for xi, yi_list, zi_list in zip(all_x, all_y, all_z):
    all_x.extend([xi] * len(yi_list))  # 把这个组的 x 扩展成与 y/z 等长
    all_y.extend(yi_list)
    all_z.extend(zi_list)

all_x = np.array(all_x)
all_y = np.array(all_y)
all_z = np.array(all_z)

# 构建统一网格
xi = np.linspace(min(all_x), max(all_x), 100)
yi = np.linspace(min(all_y), max(all_y), 100)
XI, YI = np.meshgrid(xi, yi)

# 插值到统一网格
ZI = griddata((all_x, all_y), all_z, (XI, YI), method='linear')

# 画等高线图
plt.figure(figsize=(8, 5))
contour = plt.contourf(XI, YI, ZI, levels=20, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour plot from irregular y-data')
plt.show()