import os
import argparse
import numpy as np
from shutil import copy

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--setup", type=int, default=None, metavar="INT",
                    help="Setup multirun.")
parser.add_argument("-a", "--average", action="store_true", default=False,
                    help="Average multirun.")
args = parser.parse_args()

with open("input", "r") as f:
    lines = f.readlines()

for line in lines:
    if "num_of_trajectories" in line:
        ntraj = int(line.split()[-1])
    if "num_of_time-steps" in line:
        tsteps = int(line.split()[-1])
    if "field_modes" in line:
        F = int(line.split()[-1])

if args.setup is not None:
    try:
        os.mkdir(f"run")
    except FileExistsError:
        pass

    for i in range(args.setup):
        try:
            os.mkdir(f"run/run{i+1}")
        except FileExistsError:
            pass
        copy("input", f"run/run{i+1}/input")
        copy("fmap.x", f"run/run{i+1}/fmap.x")

elif args.average:
    with open("input", "r") as f:
        lines = f.readlines()

    for line in lines:
        if "num_of_trajectories" in line:
            ntraj = int(line.split()[-1])
        if "num_of_time-steps" in line:
            tsteps = int(line.split()[-1])
        if "field_modes" in line:
            F = int(line.split()[-1])

    dirs = []
    for root, dir, files in os.walk("run"):
        dirs += dir

    fnames = [i for i in os.listdir("run/run1") if ".out" in i]

    print(f"F {F}")
    print(f"ntraj {ntraj}")
    print(f"tsteps {tsteps}")
    print(f"file names {fnames}")
    print(f"directories {dirs}")

    dat = {}
    for f in fnames:
        s = np.genfromtxt(f"run/run1/{f}").shape
        dat[f] = np.zeros(s, dtype=np.float64)
        for r in dirs:
            dat[f] += np.genfromtxt(f"run/{r}/{f}")
        dat[f] /= float(len(dirs))
        np.savetxt(f, dat[f])
