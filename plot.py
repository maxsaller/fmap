import numpy as np
from scipy import integrate as integrate
from matplotlib import pyplot as plt

plt.rcParams.update({"text.usetex": True})

# Read data files
exact = np.genfromtxt("nph_exact")

LSC = np.genfromtxt("Npop.out")
unity = np.genfromtxt("Nuni.out")
mLSCz = np.genfromtxt("Ntmz.out")
mLSCt = np.genfromtxt("Ntmt.out")
time = LSC[:, 0]

Ndot = np.genfromtxt("Ndot.out")
NQI = np.genfromtxt("NQI.out")

# Integrate time derivative
Nint = np.zeros_like(Ndot)
for i in range(1, len(time)):
    for j in range(1, 3):
        Nint[i, j] = integrate.simps(Ndot[:i, j], time[:i])

fig, ax = plt.subplots()
ax.plot(time, Ndot[:, 1])
ax.plot(time, Ndot[:, 2])
ax.plot(time, Nint[:, 1], "C0--", label=r"$[N \otimes |1\rangle\langle 1|](t)$")
ax.plot(time, Nint[:, 2], "C1--", label=r"$[N \otimes |2\rangle\langle 2|](t)$")
ax.plot(time, Nint[:, 2] - Nint[:, 1], label=r"$[N \otimes (|2\rangle\langle 2|-|1\rangle\langle 1|)](t)$")
ax.legend(frameon=False)
plt.tight_layout()
plt.show()


Nder = np.empty_like(NQI)
Nder[:, 0] = time
for i in range(1, 3):
    Nder[:, i] = (NQI[:, i] + np.sum(Nint, axis=1)) / 2.
    # Nder[:, i] = (NQI[:, i] - Nint[:, 1] + Nint[:, 2]) / 2.
print(np.sum(Nint, axis=1))

fig, ax = plt.subplots(figsize=(12, 4), dpi=100)

ax.plot(exact[:, 0], exact[:, 1], "k", label=r"Exact")
ax.plot(time, LSC[:, 2], label=r"LSC I")
ax.plot(time, unity[:, 2], label=r"$\hat{I}\mapsto 2\phi$")
# ax.plot(time, mLSCz[:, 2], label=r"mLSC/LSC-z")
# ax.plot(time, mLSCt[:, 2], label=r"mLSC/LSC-t")
ax.plot(time, Nder[:, 2], label=r"$C_{QIN}(t)$")

ax.legend(fontsize="x-large", frameon=False, ncol=6)
ax.set_xlabel(r"time (a.u.)", fontsize="x-large")
ax.set_ylabel(r"Average photon number", fontsize="x-large")

plt.tight_layout()
plt.savefig("N_CQIN.png", format="png", dpi=300)
plt.show()
