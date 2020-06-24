import numpy as np
# from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# exact = np.genfromtxt("exact")
# Cpop = np.genfromtxt("Cpop.out")
# Cimp = np.genfromtxt("Cimp.out")
# fig, ax = plt.subplots()
# ax.plot(exact[:, 0], exact[:, 1], "k")
# ax.plot(exact[:, 0], exact[:, 2], "k--")
# ax.plot(Cpop[:, 0], Cpop[:, 3], "C0")
# ax.plot(Cpop[:, 0], Cpop[:, 4], "C0--")
# ax.plot(Cimp[:, 0], Cimp[:, 3], "C1")
# ax.plot(Cimp[:, 0], Cimp[:, 4], "C1--")
# plt.tight_layout()
# plt.show(block=False)

omega = np.genfromtxt("freq.out")
I_pop = np.genfromtxt("I_pop.out")
NP_pop = np.genfromtxt("Nph_pop.out")

time = NP_pop[:, 0]
SNP_pop = NP_pop[:, 1]
NP_pop_m = NP_pop[:, 2:]
dt = time[1] - time[0]

time2 = I_pop[:, 0]
L = 236215.76557822127
r = np.arange(101) * 0.01 * L

# print(f"time: {len(time)}")
# print(f"time2: {len(time2)}")

fig, ax = plt.subplots(3, 1, figsize=(8, 12), dpi=100)
ln0, = ax[0].plot(omega, NP_pop_m[0])
ln1, = ax[1].plot(time, SNP_pop, "C0")
ln2, = ax[2].plot(r, I_pop[0, 1:])
mrk0, = ax[0].plot([0.394, 0.394],
                   [1.1 * np.min(NP_pop_m), 1.05 * np.max(NP_pop_m)],
                   'red', lw=2)
mrk1, = ax[1].plot([0, 0], [0, np.max(SNP_pop)], 'red', lw=2)


def init():
    ax[0].set_xlim([0, np.max(omega)])
    ax[0].set_ylim([1.1 * np.min(NP_pop_m), 1.1 * np.max(NP_pop_m)])
    ax[0].set_ylabel(r"Number of photons")
    ax[0].set_xlabel(r"$\omega\,/\,a.u.$")
    ax[1].set_xlim([0, np.max(time)])
    ax[1].set_ylim([-0.1, 1.1 * np.max([np.max(SNP_pop), np.max(SNP_pop)])])
    ax[1].set_xlabel(r"$t\,/\,a.u.$")
    ax[1].set_ylabel(r"Number of photons")
    ax[2].set_xlim([0, L])
    ax[2].set_ylim([-1e-7, 1e-5])
    ax[2].set_xlabel(r"$r\,/\,a.u.$")
    ax[2].set_ylabel(r"Intensity / a.u.")
    plt.subplots_adjust(right=0.97, bottom=0.05, top=1.0)
    return ln0, mrk0, ln1, mrk1, ln2


def animate(frame):
    ln0.set_data(omega, NP_pop_m[frame])
    mrk0.set_data([0.394, 0.394],
                  [1.1 * np.min(NP_pop_m), 1.05 * np.max(NP_pop_m)])
    mrk1.set_data([frame * dt, frame * dt], [0, np.max(SNP_pop)])
    if frame % ((len(time) - 1) / 100) == 0:
        print(frame)
        ln2.set_data(r, I_pop[frame // 210, 1:])
    return ln0, mrk0, ln1, mrk1, ln2


A = ani.FuncAnimation(fig, animate, frames=range(0, len(time), 21),
                      init_func=init, blit=True, interval=50. / len(time),
                      repeat=True, repeat_delay=5000)

plt.show(block=True)
