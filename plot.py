import matplotlib
import numpy as np
# from polylib import units
import matplotlib.pyplot as plt
# import matplotlib.animation as ani

matplotlib.use("TkAgg")

raw1 = np.genfromtxt("EFI_init1.out")
raw2 = np.genfromtxt("EFI_init2.out")

time = raw1[:, 0]
EFI1 = raw1[:, 1:]
EFI2 = raw2[:, 1:]
r = 12.5 / EFI1.shape[-1] * np.arange(EFI1.shape[-1], dtype=np.float64)

fig, ax = plt.subplots(EFI1.shape[0], figsize=(12, 12))

for i in range(EFI1.shape[0]):
    ax[i].plot(r, EFI1[i])
    ax[i].plot(r, EFI2[i])
    # ax[i].set_ylim([-1e-4, 6e-3])
    ax[i].set_xticks([])

plt.tight_layout(h_pad=-0.3, w_pad=1)
plt.show()


# fig, ax = plt.subplots()
# ax.plot(exact[:, 0], exact[:, 1], "k")
# ax.plot(exact[:, 0], exact[:, 2], "k--")
# ax.plot(Cpop[:, 0], Cpop[:, 3], "C0")
# ax.plot(Cpop[:, 0], Cpop[:, 4], "C0--")
# ax.plot(Cimp[:, 0], Cimp[:, 3], "C1")
# ax.plot(Cimp[:, 0], Cimp[:, 4], "C1--")
# plt.tight_layout()
# plt.show(block=False)

# omega = np.genfromtxt("freq.out")
# I_pop = np.genfromtxt("I_pop.out")
# NP_pop = np.genfromtxt("Nph_pop.out")
# exact = np.genfromtxt("exact")
# Cpop = np.genfromtxt("Cpop.out")
# Cimp = np.genfromtxt("Cimp.out")

# time = NP_pop[:, 0]
# SNP_pop = NP_pop[:, 1]
# NP_pop_m = NP_pop[:, 2:]
# dt = time[1] - time[0]

# time2 = I_pop[:, 0]
# L = 236215.76557822127
# r = units.Length(np.arange(I_pop.shape[-1] - 1) * 0.001 * L, "au").get("m")
# r *= 1e6

# fig, ax = plt.subplots(3, 1, figsize=(22, 10), dpi=100)
# ax[0].bar(omega, NP_pop_m[0], width=0.9 * (omega[1] - omega[0]))
# ln1, = ax[1].plot(time, SNP_pop, "C0")
# ln2, = ax[2].plot(r, I_pop[0, 1:])
# mrk0, = ax[0].plot([0.394, 0.394],
#                    [1.1 * np.min(NP_pop_m), 1.05 * np.max(NP_pop_m)],
#                    'red', lw=2)
# mrk1, = ax[1].plot([0, 0], [0, np.max(SNP_pop)], 'red', lw=2)
# ax[1].plot(time, Cpop[:, 3], "C1--")
# ax[1].plot(time, Cpop[:, 4], "C1")
# # ax[1].plot(time, Cimp[:, 3], "C2--")
# # ax[1].plot(time, Cimp[:, 4], "C2")


# def init():
#     ax[0].set_xlim([0, np.max(omega)])
#     ax[0].set_ylim([1.1 * np.min(NP_pop_m), 1.1 * np.max(NP_pop_m)])
#     ax[0].set_ylabel(r"Number of photons")
#     ax[0].set_xlabel(r"$\omega\,/\,a.u.$")
#     ax[1].set_xlim([0, np.max(time)])
#     # ax[1].set_ylim([-0.1, 1])
#     ax[1].set_xlabel(r"$t\,/\,a.u.$")
#     ax[1].set_ylabel(r"Number of photons")
#     ax[2].set_xlim([0, r[-1]])
#     ax[2].set_ylim([-1e-4, 1e-3])
#     ax[2].set_xlabel(r"$r\,/\,\mu m$")
#     ax[2].set_ylabel(r"Intensity / a.u.")
#     plt.subplots_adjust(left=0.05, right=0.99, bottom=0.05, top=0.99)


# def animate(frame):
#     for bar in ax[0].containers:
#         bar.remove()
#     ax[0].bar(omega, NP_pop_m[frame], width=0.9 * (omega[1] - omega[0]),
#               color="C0")
#     mrk0.set_data([0.394, 0.394],
#                   [1.1 * np.min(NP_pop_m), 1.05 * np.max(NP_pop_m)])
#     mrk1.set_data([frame * dt, frame * dt], [0, np.max(SNP_pop)])
#     if frame % ((len(time) - 1) / 100) == 0:
#         ln2.set_data(r, I_pop[frame // 210, 1:])


# A = ani.FuncAnimation(fig, animate, frames=range(0, len(time), 21),
#                       interval=(1000 / 24), init_func=init)
# print("Saving animation!")
# A.save("ani.mp4", writer="ffmpeg")
# print("Done!")
# plt.show(block=True)
