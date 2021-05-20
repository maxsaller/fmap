import numpy as np
from matplotlib import pyplot as plt

N = 10

shp = np.shape(np.genfromtxt("run/run1/Npop.out"))
Npopraw = np.empty((N, *shp), dtype=np.float64)
for i in range(N):
    Npopraw[i] = np.genfromtxt(f"run/run{i+1}/Npop.out")
Npop = np.mean(Npopraw, axis=0)
Npopstd = np.std(Npopraw, axis=0)

NQIraw = np.empty((N, *shp), dtype=np.float64)
for i in range(N):
    NQIraw[i] = np.genfromtxt(f"run/run{i+1}/NQI.out")
NQI = np.mean(NQIraw, axis=0)
NQIstd = np.std(NQIraw, axis=0)

time = Npopraw[0, :, 0]

exact = np.genfromtxt("nph_exact")

fig, ax = plt.subplots(2, sharex=True, sharey=False)

ax[0].plot(exact[:, 0], exact[:, 1], "black")
ax[0].fill_between(time, y1=Npop[:, 3] + Npopstd[:, 3], y2=Npop[:, 3] - Npopstd[:, 3],
                   color="lightgray")
ax[0].plot(time, Npop[:, 3], "C0")

ax[1].plot(exact[:, 0], exact[:, 1], "black")
ax[1].fill_between(time, y1=NQI[:, 3] + NQIstd[:, 3], y2=NQI[:, 3] - NQIstd[:, 3],
                   color="lightgray")
ax[1].plot(time, NQI[:, 3], "C1")

ax[0].set_ylabel("Photon number")
ax[1].set_ylabel("Photon number")
ax[1].set_xlabel("time (a.u.)")

plt.tight_layout()
plt.show()
