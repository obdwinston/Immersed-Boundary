import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Polygon

# =============================================================================
# CONFIGURATION
# =============================================================================

DATA_DIR = os.path.join("outputs", "data")
FORCE_FILE = os.path.join(DATA_DIR, "forces.dat")
BODY_FILE = "body.dat"
ANIMATION_FILE = os.path.join("outputs", "animation.mp4")
METRICS_FILE = os.path.join("outputs", "metrics.dat")

Lc = 1.0
U_inf = 1.0
x_view = (-2, 6)
y_view = (-1.5, 1.5)

# =============================================================================
# DATA LOADING
# =============================================================================


def load_body(filename):
    """Load body geometry and normalise to unit chord centred at origin."""
    raw = np.loadtxt(filename)

    # Compute arc-length weighted centroid
    arc = np.hypot(
        np.diff(raw[:, 0], append=raw[0, 0]),
        np.diff(raw[:, 1], append=raw[0, 1]),
    )
    arc_total = np.sum(arc)
    xc = np.sum(raw[:, 0] * arc) / arc_total
    yc = np.sum(raw[:, 1] * arc) / arc_total

    # Centre and normalise
    local = raw.copy()
    local[:, 0] -= xc
    local[:, 1] -= yc
    dists = np.hypot(
        local[:, 0, None] - local[:, 0],
        local[:, 1, None] - local[:, 1],
    )
    local /= np.max(dists)

    return local


def load_frames(data_dir):
    """Load frame metadata and grid from data files."""
    files = sorted(glob.glob(os.path.join(data_dir, "fields_*.dat")))
    if not files:
        raise FileNotFoundError(f"No files found in {data_dir}/")

    frames = []
    for f in files:
        with open(f) as fh:
            fh.readline()
            header = fh.readline().split()
            nx, ny = int(header[0]), int(header[1])
        frames.append(
            dict(
                file=f,
                t=float(header[2]),
                Xc=float(header[3]),
                Yc=float(header[4]),
                Th=float(header[5]),
            )
        )

    # Load grid from last file
    data = np.loadtxt(files[-1], comments="#", skiprows=2)
    x = data[:, 0].reshape(ny, nx, order="C")
    y = data[:, 1].reshape(ny, nx, order="C")
    dx = x[0, 1] - x[0, 0]
    dy = y[1, 0] - y[0, 0]

    return frames, x, y, nx, ny, dx, dy


def load_field(filename, nx, ny):
    """Load velocity field and compute magnitude."""
    data = np.loadtxt(filename, comments="#", skiprows=2)
    u = data[:, 2].reshape(ny, nx, order="C")
    v = data[:, 3].reshape(ny, nx, order="C")
    return np.hypot(u, v), u, v


def load_forces(filename):
    """Load force coefficients from file."""
    if os.path.exists(filename):
        return np.loadtxt(filename, comments="#")
    return None


# =============================================================================
# METRICS
# =============================================================================


def compute_metrics(forces):
    """Compute force coefficient statistics from second half of data."""
    if forces is None or len(forces) <= 20:
        return None

    # Use second half to skip transient
    n = len(forces) // 2
    d = forces[n:]
    t, cd, cl, cm = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    dt = t[1] - t[0]

    # Find zero crossings for analysis window
    crossings = np.where((cl[:-1] <= 0) & (cl[1:] > 0))[0]
    if len(crossings) < 2:
        return None

    start, end = crossings[0], crossings[-1]
    t_win = t[start : end + 1]
    cl_win = cl[start : end + 1]
    cd_win = cd[start : end + 1]
    cm_win = cm[start : end + 1]

    # FFT to find dominant frequency
    spec = np.abs(np.fft.rfft(cl_win - np.mean(cl_win)))
    freqs = np.fft.rfftfreq(len(cl_win), d=dt)
    peak = np.argmax(spec[1:]) + 1
    freq_dominant = freqs[peak]
    cycles = int(round((t_win[-1] - t_win[0]) * freq_dominant))

    return {
        "window": f"t = {t_win[0]:.1f} - {t_win[-1]:.1f}",
        "cycles": cycles,
        "St": freq_dominant * Lc / U_inf,
        "CL_avg": np.mean(cl_win),
        "CL_min": np.min(cl_win),
        "CL_max": np.max(cl_win),
        "CL_amp": np.max(cl_win) - np.min(cl_win),
        "CD_avg": np.mean(cd_win),
        "CD_min": np.min(cd_win),
        "CD_max": np.max(cd_win),
        "CD_amp": np.max(cd_win) - np.min(cd_win),
        "CM_avg": np.mean(cm_win),
        "CM_min": np.min(cm_win),
        "CM_max": np.max(cm_win),
        "CM_amp": np.max(cm_win) - np.min(cm_win),
    }


def print_metrics(metrics):
    """Print metrics summary to console."""
    print(f"{metrics['cycles']} complete cycles ({metrics['window']})")
    print(f"St = {metrics['St']:.4f}")
    print(
        f"CL: avg={metrics['CL_avg']:.4f}, min={metrics['CL_min']:.4f}, "
        f"max={metrics['CL_max']:.4f}, amp={metrics['CL_amp']:.4f}"
    )
    print(
        f"CD: avg={metrics['CD_avg']:.4f}, min={metrics['CD_min']:.4f}, "
        f"max={metrics['CD_max']:.4f}, amp={metrics['CD_amp']:.4f}"
    )
    print(
        f"CM: avg={metrics['CM_avg']:.4f}, min={metrics['CM_min']:.4f}, "
        f"max={metrics['CM_max']:.4f}, amp={metrics['CM_amp']:.4f}"
    )


def save_metrics(metrics, filename):
    """Save metrics to file."""
    keys = [
        "St",
        "CL_avg",
        "CL_min",
        "CL_max",
        "CL_amp",
        "CD_avg",
        "CD_min",
        "CD_max",
        "CD_amp",
        "CM_avg",
        "CM_min",
        "CM_max",
        "CM_amp",
    ]
    with open(filename, "w") as f:
        f.write(f"# Analysis window: {metrics['window']}\n")
        f.write(f"# Complete cycles: {metrics['cycles']}\n")
        for key in keys:
            f.write(f"{key} {metrics[key]:.4f}\n")


# =============================================================================
# UTILITIES
# =============================================================================


def transform_body(body_local, Xc, Yc, Th):
    """Transform body from local to world coordinates."""
    cos_th, sin_th = np.cos(Th), np.sin(Th)
    x_world = Xc + Lc * (body_local[:, 0] * cos_th - body_local[:, 1] * sin_th)
    y_world = Yc + Lc * (body_local[:, 0] * sin_th + body_local[:, 1] * cos_th)
    return np.column_stack([x_world, y_world])


def compute_colour_limits(frames, nx, ny, dx, dy):
    """Compute colour limits from second half of simulation."""
    t_mid = frames[-1]["t"] / 2
    i_mid = next(i for i, f in enumerate(frames) if f["t"] >= t_mid)

    vmax_u, vmax_w = 0.0, 0.0
    for i in range(i_mid, len(frames)):
        U_i, u, v = load_field(frames[i]["file"], nx, ny)
        w = np.gradient(v, dx, axis=1) - np.gradient(u, dy, axis=0)
        vmax_u = max(vmax_u, np.nanmax(U_i))
        vmax_w = max(vmax_w, np.nanmax(np.abs(w)))

    return vmax_u, vmax_w * 0.4


# =============================================================================
# ANIMATION
# =============================================================================


def setup_figure(frames, forces, x, y, nx, ny, dx, dy, body_local, vmax_u, vmax_w):
    """Set up the animation figure with three subplots."""
    aspect_ratio = (x_view[1] - x_view[0]) / (y_view[1] - y_view[0])
    fig, (ax_u, ax_w, ax_force) = plt.subplots(3, 1, figsize=(10, 10))

    # Initial frame data
    f0 = frames[0]
    U0, u0, v0 = load_field(f0["file"], nx, ny)
    w0 = np.gradient(v0, dx, axis=1) - np.gradient(u0, dy, axis=0)
    body_xy = transform_body(body_local, f0["Xc"], f0["Yc"], f0["Th"])

    # Velocity magnitude plot
    pcm_u = ax_u.pcolormesh(
        x, y, U0, shading="auto", cmap="RdYlBu_r", vmin=0, vmax=vmax_u
    )
    cax_u = ax_u.inset_axes([1.02, 0, 0.02, 1])
    fig.colorbar(pcm_u, cax=cax_u, label="|velocity|")
    title = ax_u.set_title(f"t = {f0['t']:.3f}", fontsize=12)
    patch_u = Polygon(body_xy, closed=True, fc="grey", ec="black", lw=1.5, zorder=10)
    ax_u.add_patch(patch_u)
    ax_u.set(xlabel="x", ylabel="y", xlim=x_view, ylim=y_view, aspect="equal")

    # Vorticity plot
    pcm_w = ax_w.pcolormesh(
        x, y, w0, shading="auto", cmap="RdBu_r", vmin=-vmax_w, vmax=vmax_w
    )
    cax_w = ax_w.inset_axes([1.02, 0, 0.02, 1])
    fig.colorbar(pcm_w, cax=cax_w, label="vorticity")
    patch_w = Polygon(body_xy, closed=True, fc="grey", ec="black", lw=1.5, zorder=10)
    ax_w.add_patch(patch_w)
    ax_w.set(xlabel="x", ylabel="y", xlim=x_view, ylim=y_view, aspect="equal")

    # Force coefficient plot
    time_marker = ax_force.axvline(0, color="k", ls="--", lw=1)
    ax_force.set(
        xlabel="Time", ylabel="Force coefficients", title="$C_D$, $C_L$, $C_M$"
    )
    if forces is not None and len(forces) > 0:
        ax_force.plot(forces[:, 0], forces[:, 1], "b-", lw=1, label="$C_D$")
        ax_force.plot(forces[:, 0], forces[:, 2], "r-", lw=1, label="$C_L$")
        ax_force.plot(forces[:, 0], forces[:, 3], "g-", lw=1, label="$C_M$")
        ax_force.grid(True, alpha=0.3)
        ax_force.legend(loc="upper right", fontsize=9)
        # Set y-limits from second half
        n_half = len(forces) // 2
        f_half = forces[n_half:]
        f_min = min(np.min(f_half[:, 1]), np.min(f_half[:, 2]), np.min(f_half[:, 3]))
        f_max = max(np.max(f_half[:, 1]), np.max(f_half[:, 2]), np.max(f_half[:, 3]))
        margin = (f_max - f_min) * 0.1
        ax_force.set_ylim(f_min - margin, f_max + margin)
    ax_force.set_box_aspect(1 / aspect_ratio)

    # Centre all axes
    for ax in (ax_u, ax_w, ax_force):
        ax.set_anchor("C")
    fig.tight_layout()

    return fig, pcm_u, pcm_w, title, time_marker, patch_u, patch_w


def create_update_function(
    frames,
    nx,
    ny,
    dx,
    dy,
    body_local,
    pcm_u,
    pcm_w,
    title,
    time_marker,
    patch_u,
    patch_w,
):
    """Create the animation update function."""
    n_frames = len(frames)

    def update_frame(frame):
        print(f"\rRendering frame {frame + 1}/{n_frames}...", end="", flush=True)
        f = frames[frame]
        U_i, u, v = load_field(f["file"], nx, ny)
        w_i = np.gradient(v, dx, axis=1) - np.gradient(u, dy, axis=0)

        pcm_u.set_array(U_i.ravel())
        pcm_w.set_array(w_i.ravel())
        title.set_text(f"t = {f['t']:.3f}")
        time_marker.set_xdata([f["t"], f["t"]])

        body_xy = transform_body(body_local, f["Xc"], f["Yc"], f["Th"])
        patch_u.set_xy(body_xy)
        patch_w.set_xy(body_xy)

        return [pcm_u, pcm_w, title, time_marker, patch_u, patch_w]

    return update_frame


# =============================================================================
# MAIN
# =============================================================================


def main():
    # Load data
    body_local = load_body(BODY_FILE)
    frames, x, y, nx, ny, dx, dy = load_frames(DATA_DIR)
    forces = load_forces(FORCE_FILE)
    print(f"{len(frames)} frames loaded (nx = {nx}, ny = {ny})")

    # Compute and print metrics
    metrics = compute_metrics(forces)
    if metrics:
        print_metrics(metrics)
        save_metrics(metrics, METRICS_FILE)
        print(f"Metrics saved to {METRICS_FILE}")

    # Compute colour limits
    print("Computing colour limits...", end="", flush=True)
    vmax_u, vmax_w = compute_colour_limits(frames, nx, ny, dx, dy)
    print(" done!")

    # Setup figure
    fig, pcm_u, pcm_w, title, time_marker, patch_u, patch_w = setup_figure(
        frames, forces, x, y, nx, ny, dx, dy, body_local, vmax_u, vmax_w
    )

    # Create animation
    update_frame = create_update_function(
        frames,
        nx,
        ny,
        dx,
        dy,
        body_local,
        pcm_u,
        pcm_w,
        title,
        time_marker,
        patch_u,
        patch_w,
    )
    anim = animation.FuncAnimation(
        fig, update_frame, frames=len(frames), interval=150, blit=False
    )
    anim.save(ANIMATION_FILE, writer="ffmpeg", fps=8, dpi=120)
    print(" done!")
    print(f"Animation saved to {ANIMATION_FILE}")
    plt.close(fig)


if __name__ == "__main__":
    main()
