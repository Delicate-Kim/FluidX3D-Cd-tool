# drop-in replacement
import argparse, glob, os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from itertools import count  # <-- keep update() running

def find_cd(path_hint):
    if path_hint and os.path.isfile(path_hint):
        return path_hint
    base = path_hint if path_hint and os.path.isdir(path_hint) else os.getcwd()
    candidates = glob.glob(os.path.join(base, "FP*/**/Cd.dat"), recursive=True)
    if not candidates:
        raise FileNotFoundError("Could not find Cd.dat. Pass --file <path/to/Cd.dat> or run from the output dir.")
    return max(candidates, key=os.path.getmtime)

def read_cd(cd_path):
    t, cd = [], []
    with open(cd_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    t.append(float(parts[0])); cd.append(float(parts[1]))
                except ValueError:
                    pass
    return t, cd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", "-f", default="", help="Path to Cd.dat (or an output directory that contains it).")
    ap.add_argument("--interval", "-i", type=int, default=500, help="Refresh interval (ms).")
    args = ap.parse_args()

    cd_path = find_cd(args.file if args.file else os.getcwd())
    print(cd_path)

    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=1.5)

    COLOR = "#14569c"
    last_marker, = ax.plot(
        [], [], linestyle='None', marker='o', ms=6,
        markerfacecolor=COLOR, markeredgecolor=COLOR, markeredgewidth=1.8,
        zorder=5
    )

    # --- Label to the LEFT of the marker ---
    last_label = ax.annotate(
        "", xy=(0, 0), xytext=(-10, 8), textcoords="offset points",
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc="w", ec="0.6", alpha=0.9),
        zorder=6
    )

    ax.set_xlabel("t (LBM timesteps)")
    ax.set_ylabel("C_d")
    ax.grid(True)
    ax.margins(x=0.02, y=0.12)
    try:
        fig.canvas.manager.set_window_title(f"Live Cd — {cd_path}")
    except Exception:
        pass

    def init():
        line.set_data([], [])
        last_marker.set_data([], [])
        last_label.set_text("")
        return line, last_marker, last_label

    def update(_frame):
        t, cd = read_cd(cd_path)
        if not t:
            return line, last_marker, last_label

        line.set_data(t, cd)
        ax.relim(); ax.autoscale_view()

        x, y = t[-1], cd[-1]
        last_marker.set_data([x], [y])

        # keep label left of the point; nudge vertically so it doesn't overlap
        ymin, ymax = ax.get_ylim()
        dy = 0.02 * (ymax - ymin)
        last_label.xy = (x, y)
        last_label.set_position((-10, 8))   # (-x, +y) offset in points (left & slightly up)
        last_label.set_text(f"C_d={y:.3f}")

        # (optional) mirror in window title
        try:
            fig.canvas.manager.set_window_title(f"Live Cd — {cd_path}  |  last C_d={y:.3f} @ t={x:g}")
        except Exception:
            pass

        return line, last_marker, last_label

    # Keep a strong ref; provide frames; disable frame caching
    ani = FuncAnimation(fig, update, init_func=init, frames=count(),
                        interval=args.interval, blit=False, cache_frame_data=False)
    fig._ani = ani

    plt.show()

if __name__ == "__main__":
    main()
