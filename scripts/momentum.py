import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 1. データの読み込み
data_file = "../data/momentum_density.csv"
try:
    density_data = np.loadtxt(data_file, delimiter=",")
    density_data = np.nan_to_num(density_data, nan=0.0, posinf=0.0, neginf=0.0)
except Exception as e:
    print(f"読み込み失敗: {e}")
    exit()

time_steps, spatial_points = density_data.shape

# 2. 運動量(p)軸の設定 (Rust側と計算式を合わせる)
N = spatial_points
dx = 0.1
dp = 2.0 * np.pi / (N * dx)
# -N/2 から N/2-1 までの配列を作り dp を掛ける
p = dp * (np.arange(N) - N / 2.0)

# 3. グラフの初期設定
fig, ax = plt.subplots(figsize=(8, 5))

# 位置空間とは違い、横軸は p になる
line, = ax.plot(p, density_data[0], color='red', lw=2, label=r"$|\phi(p)|^2$")

max_val = np.max(density_data)
if max_val <= 0: max_val = 1.0

# 運動量の広がりを見るため、X軸の範囲を適度に絞る
ax.set_xlim(-10, 10) 
ax.set_ylim(0, max_val * 1.1)
ax.set_title("Time Evolution of Momentum Distribution")
ax.set_xlabel("Momentum (p)")
ax.set_ylabel("Probability Density")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.6)

# 4. アニメーション更新
def update(frame):
    line.set_ydata(density_data[frame])
    return line,

ani = animation.FuncAnimation(fig, update, frames=time_steps, interval=30, blit=True)
plt.show()