import numpy as np
from matplotlib import pyplot as plt

plt.rcParams.update({"text.usetex": True})

dat = np.genfromtxt("Ipop.out")

L = 12.5
x = np.arange(1000) * L / 1000

# for i in range(22):
#     plt.plot(x, dat[i, 1:])
#     plt.ylim([-1e-4, 6e-4])
#     plt.show()

fig, ax = plt.subplots(6, figsize=(8, 6), sharex=True, sharey=True)

phax = fig.add_subplot(111, frameon=False)
phax.tick_params(labelcolor="none", top=False, bottom=False, left=False,
                 right=False)

ax[0].plot(x, dat[1, 1001:])
ax[1].plot(x, dat[6, 1001:])
ax[2].plot(x, dat[11, 1001:])
ax[3].plot(x, dat[17, 1001:])
ax[4].plot(x, dat[18, 1001:])
ax[5].plot(x, dat[21, 1001:])


ax[0].set_xlim(0, 12.5)
# ax[0].set_ylim(-1e-3, 5e-3)
for i in range(6):
    ax[i].tick_params(axis="both", which="both", direction="in",
                      bottom=True, top=True, left=True, right=True)
    ax[i].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    if i != 0:
        ax[i].yaxis.offsetText.set_visible(False)

ax[0].text(0.25, 0.0035, r"$t=100\,\mathrm{a.u.}$")
ax[1].text(0.25, 0.0035, r"$t=600\,\mathrm{a.u.}$")
ax[2].text(0.25, 0.0035, r"$t=1100\,\mathrm{a.u.}$")
ax[3].text(0.25, 0.0035, r"$t=1700\,\mathrm{a.u.}$")
ax[4].text(0.25, 0.0035, r"$t=1800\,\mathrm{a.u.}$")
ax[5].text(0.25, 0.0035, r"$t=2100\,\mathrm{a.u.}$")

ax[5].set_xlabel(r"Cavity position (\textmu m)", fontsize="large")
phax.set_ylabel(r"Intensity", fontsize="large")

plt.tight_layout(h_pad=0.25, w_pad=0)
# plt.savefig("2D_int_ehr.eps", formar="eps", dpi=300)
plt.show()
