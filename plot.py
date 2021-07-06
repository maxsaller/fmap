import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from scipy import integrate

# N = 10

# shp = np.shape(np.genfromtxt("run/run1/Npop.out"))
# Npopraw = np.empty((N, *shp), dtype=np.float64)
# for i in range(N):
#     Npopraw[i] = np.genfromtxt(f"run/run{i+1}/Npop.out")
# Npop = np.mean(Npopraw, axis=0)
# Npopstd = np.std(Npopraw, axis=0)

# Npop2raw = np.empty((N, *shp), dtype=np.float64)
# for i in range(N):
#     Npop2raw[i] = np.genfromtxt(f"run2/run{i+1}/Npop.out")
# Npop2 = np.mean(Npop2raw, axis=0)
# Npop2std = np.std(Npop2raw, axis=0)

# NQIraw = np.empty((N, *shp), dtype=np.float64)
# for i in range(N):
#     NQIraw[i] = np.genfromtxt(f"run/run{i+1}/NQI.out")
# NQI = np.mean(NQIraw, axis=0)
# NQIstd = np.std(NQIraw, axis=0)

# time = Npopraw[0, :, 0]

# exact = np.genfromtxt("nph_exact")

# fig, ax = plt.subplots(3, sharex=True, sharey=False)

# ax[0].plot(exact[:, 0], exact[:, 1], "black")
# ax[0].fill_between(time, y1=Npop[:, 3] + Npopstd[:, 3], y2=Npop[:, 3] - Npopstd[:, 3],
#                    color="lightgray")
# ax[0].plot(time, Npop[:, 3], "C1")

# ax[1].plot(exact[:, 0], exact[:, 1], "black")
# ax[1].fill_between(time, y1=Npop2[:, 3] + Npop2std[:, 3], y2=Npop2[:, 3] - Npop2std[:, 3],
#                    color="lightgray")
# ax[1].plot(time, Npop2[:, 3], "C2")

# ax[2].plot(exact[:, 0], exact[:, 1], "black")
# ax[2].fill_between(time, y1=NQI[:, 3] + NQIstd[:, 3], y2=NQI[:, 3] - NQIstd[:, 3],
#                    color="lightgray")
# ax[2].plot(time, NQI[:, 3], "C3")

# ax[0].set_ylabel("Photon number")
# ax[1].set_ylabel("Photon number")
# ax[2].set_ylabel("Photon number")
# ax[2].set_xlabel("time (a.u.)")

# plt.tight_layout()
# plt.show()

# N = 10

# shp = np.shape(np.genfromtxt("run/run1/Cpop.out"))
# Cpopraw = np.empty((N, *shp), dtype=np.float64)
# for i in range(N):
#     Cpopraw[i] = np.genfromtxt(f"run/run{i+1}/Cpop.out")
# Cpop = np.mean(Cpopraw, axis=0)
# Cpopstd = np.std(Cpopraw, axis=0)

# Cpop2raw = np.empty((N, *shp), dtype=np.float64)
# for i in range(N):
#     Cpop2raw[i] = np.genfromtxt(f"run2/run{i+1}/Cpop.out")
# Cpop2 = np.mean(Cpop2raw, axis=0)
# Cpop2std = np.std(Cpop2raw, axis=0)

# Cimpraw = np.empty((N, *shp), dtype=np.float64)
# for i in range(N):
#     Cimpraw[i] = np.genfromtxt(f"run/run{i+1}/Cimp.out")
# Cimp = np.mean(Cimpraw, axis=0)
# Cimpstd = np.std(Cimpraw, axis=0)


# exact11 = np.genfromtxt("exact11")
# exact22 = np.genfromtxt("exact22")
# exact33 = np.genfromtxt("exact33")

# fig, ax = plt.subplots(3, sharex=True, sharey=False)

# for a in ax:
#     a.plot(exact11[:, 0], exact11[:, 1], "k")
#     a.plot(exact22[:, 0], exact22[:, 1], "k--")
#     a.plot(exact33[:, 0], exact33[:, 1], "k:")

# ax[0].plot(time, Cpop[:, 7], "C1")
# ax[0].plot(time, Cpop[:, 8], "C1--")
# ax[0].plot(time, Cpop[:, 9], "C1:")

# ax[1].plot(time, Cpop2[:, 7], "C1")
# ax[1].plot(time, Cpop2[:, 8], "C1--")
# ax[1].plot(time, Cpop2[:, 9], "C1:")

# ax[2].plot(time, Cimp[:, 7], "C3")
# ax[2].plot(time, Cimp[:, 8], "C3--")
# ax[2].plot(time, Cimp[:, 9], "C3:")

# ax[0].set_ylabel("Atomic population")
# ax[1].set_ylabel("Atomic population")
# ax[2].set_ylabel("Atomic population")
# ax[2].set_xlabel("time (a.u.)")

# plt.tight_layout()
# plt.show()

cav_steps = 1000
L = 236215.76557822127
r = np.arange(cav_steps) * L / float(cav_steps - 1)
L = 12.5
rmum = np.arange(cav_steps) * L / float(cav_steps - 1)

shape = np.genfromtxt("run/run1/Ipop.out").shape
raw = np.empty((5, *shape), dtype=np.float64)
for i in range(5):
    raw[i] = np.genfromtxt(f"run/run{i + 1}/Ipop.out")
Int = np.mean(raw, axis=0)

fig, ax = plt.subplots()
ln1, = ax.plot(rmum, 3.e-3 * np.ones_like(r), "C3--")
ln, = ax.plot(rmum, Int[0, 2001:], "C0")


def init():
    ax.set_xlim([0, rmum[-1]])
    ax.set_ylim([-1e-3, 1e-2])
    ax.set_xlabel(r"$r\,(\mu m)$", fontsize="large")
    ax.set_ylabel(r"$I(r)$", fontsize="large")
    return ln, ln1,


def animate(ts):
    ln.set_data(rmum, Int[ts, 2001:])
    return ln, ln1,


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=211, blit=True, interval=50)
plt.show()

fig, ax = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

exact = np.genfromtxt("nph_exact")

N = np.empty(Int.shape[0], dtype=np.float64)
for i in range(Int.shape[0]):
    N[i] = integrate.trapezoid(Int[i, 2001:], x=r)

ax[0].plot(Int[:, 0], N)
ax[0].set_ylim([-1, 19])

ax[1].plot(exact[:, 0], exact[:, 1], "k", label="Exact")
ax[1].plot(Int[:, 0], N / 8., "C1", label=r"$\frac{1}{8}\int I_{\mathrm{LSC\,I}}(r) \mathrm{d}r$")
ax[1].set_ylim([-0.1, 2.5])

ax[0].set_xlabel("time (a.u.)", fontsize="large")
ax[1].set_xlabel("time (a.u.)", fontsize="large")
ax[0].set_ylabel("Photon number", fontsize="large")

ax[1].legend(frameon=False, fontsize="xx-large")

plt.tight_layout()
plt.show()
