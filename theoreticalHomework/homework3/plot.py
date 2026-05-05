import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# --- 1. 定义 B-spline 函数 ---

def B2_uniform_scalar(x, i):
    """
    计算 B_i^2(x) 在 t_j = j (均匀节点) 时的标量值。
    支撑区间为 [i-1, i+2]
    """
    # 区间 1: [i-1, i]
    if (x >= i - 1) and (x < i):
        return 0.5 * (x - i + 1)**2
    # 区间 2: [i, i+1]
    # B_i^2(x) = 0.5 * [ (x-i+1)(i+1-x) + (i+2-x)(x-i) ]
    # 展开后: -x^2 + (2i+1)x - i^2 - i + 0.5
    elif (x >= i) and (x < i + 1):
        term1 = (x - i + 1) * (i + 1 - x)
        term2 = (i + 2 - x) * (x - i)
        return 0.5 * (term1 + term2)
    # 区间 3: [i+1, i+2]
    elif (x >= i + 1) and (x <= i + 2):
        return 0.5 * (i + 2 - x)**2
    # 支撑区间之外
    else:
        return 0.0

# 将标量函数向量化，使其可以处理 numpy 数组
# 'i' 参数保持不变 (excluded)
B2_vec = np.vectorize(B2_uniform_scalar, excluded=['i'])

# --- 2. 设置字体和颜色 ---

# 尝试加载您指定的字体
try:
    font_path_title = '/usr/share/fonts/custom/FmlSong-Black.ttf'
    font_title_prop = fm.FontProperties(fname=font_path_title, size=16)
except FileNotFoundError:
    print(f"Warning: Title font not found at {font_path_title}. Using default font.")
    font_title_prop = fm.FontProperties(size=16)

try:
    font_path_labels = '/usr/share/fonts/custom/timesbd.ttf'
    font_labels_prop = fm.FontProperties(fname=font_path_labels, size=12)
except FileNotFoundError:
    print(f"Warning: Labels font not found at {font_path_labels}. Using default font.")
    font_labels_prop = fm.FontProperties(size=12)

# 您选择的颜色
colors = ['#f9f871', '#ff9671', '#00c9a7', '#0081cf', '#d65db1', '#845ec2']
# 如果颜色不够，可以循环使用
color_list = ['#f9f871', '#ff9671', '#00c9a7', '#0081cf']


# --- 3. 绘图 ---

# 定义 x 轴范围
x = np.linspace(-2, 8, 500)

# 设置图像大小
plt.figure(figsize=(10, 6))

# 绘制 i = 0 到 5 的 B-splines
i_values = range(6)
for i in i_values:
    y = B2_vec(x, i=i)
    plt.plot(x, y, 
             label=f'$B_{i}^2(x)$', 
             color=color_list[i % len(color_list)], 
             linewidth=2)

# --- 4. 自定义和保存图像 ---

# 设置标题和坐标轴标签，使用指定字体
plt.title("Plot of B_i^2(x) for uniform knots t_i = i", fontproperties=font_title_prop)
plt.xlabel("x", fontproperties=font_labels_prop)
plt.ylabel("B_i^2(x)", fontproperties=font_labels_prop)

# 设置坐标轴范围
plt.ylim(0, 0.8) # B-spline 的最大值是 0.75
plt.xlim(-2, 8)

# 添加图例和网格
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# 设置坐标轴刻度线朝内
plt.tick_params(direction='in')

# 保存图像
output_filename = "quadratic_b_spline_plot.png"
plt.savefig(output_filename, dpi=150, bbox_inches='tight')

print(f"Plot saved to {output_filename}")
