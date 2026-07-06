import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 1. データの読み込み
data_file = "../data/probability_density.csv"
try:
    density_data = np.loadtxt(data_file, delimiter=",")
except Exception as e:
    print(f"データの読み込みに失敗しました: {e}")
    exit()

# 【重要】データ内にNaNやInfが含まれていた場合、0.0に置き換えてクラッシュを防ぐ
density_data = np.nan_to_num(density_data, nan=0.0, posinf=0.0, neginf=0.0)

time_steps, spatial_points = density_data.shape

# 2. グラフの初期設定
fig, ax = plt.subplots(figsize=(8, 5))
dx = 0.1
x = np.linspace(0, spatial_points * dx, spatial_points)

# 修正1: raw文字列 (r"...") を使用して SyntaxWarning を解消
line, = ax.plot(x, density_data[0], color='blue', lw=2, label=r"$|\psi(x)|^2$")

# 修正2: 安全なY軸の最大値設定
max_val = np.max(density_data)
if max_val <= 0:
    max_val = 1.0  # データが全て0等の異常時のフォールバック

ax.set_xlim(0, spatial_points * dx)
ax.set_ylim(0, max_val * 1.1)
ax.set_title("Time Evolution of Probability Distribution")
ax.set_xlabel("Position (x)")
ax.set_ylabel("Probability Density")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.6)

# 3. アニメーション更新用の関数
def update(frame):
    line.set_ydata(density_data[frame])
    return line,

# 4. アニメーションの生成
ani = animation.FuncAnimation(
    fig, 
    update, 
    frames=time_steps, 
    interval=30, 
    blit=True
)

plt.show()