#!/usr/bin/env python3
"""
plot_results.py

Reads the data files produced by classical_double_pendulum.f90 and
quantum_double_pendulum.f90 and produces two PNG figures:

  overview.png            - classical trajectory/phase-portrait/energy,
                             quantum expectation values, norm & energy drift,
                             and a classical-vs-quantum angle comparison.
  quantum_snapshots.png   - a strip of |psi(theta1,theta2)|^2 heatmaps at
                             several times, showing the wavepacket evolve.

Requires: numpy, matplotlib  (pip install numpy matplotlib)
Run this in the same directory as the three .dat files:
    python3 plot_results.py
"""
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. Load classical trajectory
# ----------------------------------------------------------------------
cls = np.loadtxt("classical_trajectory.dat")
t_c, th1_c, th2_c, w1_c, w2_c, x1_c, y1_c, x2_c, y2_c, E_c = cls.T

# ----------------------------------------------------------------------
# 2. Load quantum observables
# ----------------------------------------------------------------------
qm = np.loadtxt("quantum_observables.dat")
t_q, e_th1_q, e_th2_q, norm_q, E_q = qm.T

# ----------------------------------------------------------------------
# 3. Load quantum density snapshots (gnuplot-style blank-line-separated blocks)
# ----------------------------------------------------------------------
def load_snapshots(path):
    snapshots = []
    cur_t = None
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("#"):
                if rows:
                    snapshots.append((cur_t, np.array(rows)))
                    rows = []
                cur_t = float(s.split("=")[1])
            elif s == "":
                continue
            else:
                rows.append([float(v) for v in s.split()])
    if rows:
        snapshots.append((cur_t, np.array(rows)))
    return snapshots

snaps = load_snapshots("quantum_density_snapshots.dat")

def snapshot_to_grid(block):
    th1 = np.unique(block[:, 0])
    th2 = np.unique(block[:, 1])
    dens = block[:, 2].reshape(len(th1), len(th2))
    return th1, th2, dens

# ----------------------------------------------------------------------
# Figure 1: overview
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# classical angles vs time
ax = axes[0, 0]
ax.plot(t_c, np.degrees(th1_c), label=r"$\theta_1$")
ax.plot(t_c, np.degrees(th2_c), label=r"$\theta_2$")
ax.set_xlabel("t (s)"); ax.set_ylabel("angle (deg)")
ax.set_title("Classical: angles vs time")
ax.legend()

# classical phase portrait
ax = axes[0, 1]
ax.plot(th1_c, w1_c, lw=0.4)
ax.set_xlabel(r"$\theta_1$ (rad)"); ax.set_ylabel(r"$\omega_1$ (rad/s)")
ax.set_title("Classical: phase portrait")

# classical bob-2 trace (the classic chaotic pattern)
ax = axes[0, 2]
ax.plot(x2_c, y2_c, lw=0.3, color="darkred")
ax.set_xlabel("x"); ax.set_ylabel("y")
ax.set_title("Classical: bob-2 trace")
ax.set_aspect("equal")

# classical energy conservation
ax = axes[1, 0]
ax.plot(t_c, E_c - E_c[0])
ax.set_xlabel("t (s)"); ax.set_ylabel(r"$E(t)-E_0$ (J)")
ax.set_title("Classical: energy conservation")

# quantum <theta1>, <theta2> vs t, overlaid with classical (different unit systems!)
ax = axes[1, 1]
ax.plot(t_q, np.degrees(e_th1_q), label=r"quantum $\langle\theta_1\rangle$")
ax.plot(t_q, np.degrees(e_th2_q), label=r"quantum $\langle\theta_2\rangle$")
ax.plot(t_c, np.degrees(th1_c), '--', alpha=0.5, label=r"classical $\theta_1$")
ax.plot(t_c, np.degrees(th2_c), '--', alpha=0.5, label=r"classical $\theta_2$")
ax.set_xlim(0, max(t_q.max(), 3))
ax.set_xlabel("t  (seconds for classical; natural units for quantum)")
ax.set_ylabel("angle (deg)")
ax.set_title("Classical vs quantum (same initial angles)\ntime axes use different unit systems")
ax.legend(fontsize=8)

# quantum norm & energy drift (diagnostics)
ax = axes[1, 2]
ax2 = ax.twinx()
l1, = ax.plot(t_q, norm_q - 1.0, color="tab:blue", label="norm drift")
l2, = ax2.plot(t_q, E_q - E_q[0], color="tab:orange", label="energy drift")
ax.set_xlabel("t (natural units)")
ax.set_ylabel("norm(t) - 1", color="tab:blue")
ax2.set_ylabel(r"$\langle H\rangle(t)-\langle H\rangle_0$", color="tab:orange")
ax.set_title("Quantum: numerical conservation check")
ax.legend(handles=[l1, l2], fontsize=8, loc="upper left")

plt.tight_layout()
plt.savefig("overview.png", dpi=150)
print("wrote overview.png")

# ----------------------------------------------------------------------
# Figure 2: quantum density snapshots
# ----------------------------------------------------------------------
n_show = min(len(snaps), 6)
idxs = np.linspace(0, len(snaps) - 1, n_show).astype(int)
fig2, axes2 = plt.subplots(1, n_show, figsize=(3 * n_show, 3.2), constrained_layout=True)
if n_show == 1:
    axes2 = [axes2]

for k, idx in enumerate(idxs):
    tt, block = snaps[idx]
    th1g, th2g, dens = snapshot_to_grid(block)
    ax = axes2[k]
    im = ax.pcolormesh(np.degrees(th1g), np.degrees(th2g), dens.T, shading="auto", cmap="inferno")
    ax.set_title(f"t = {tt:.2f}")
    ax.set_xlabel(r"$\theta_1$ (deg)")
    if k == 0:
        ax.set_ylabel(r"$\theta_2$ (deg)")

fig2.suptitle(r"Quantum density $|\psi(\theta_1,\theta_2)|^2$ over time")
plt.savefig("quantum_snapshots.png", dpi=150)
print("wrote quantum_snapshots.png")
