import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("--line","-l", type=float, default=-999999.0, help="reference line")
ap.add_argument("--shift", "-s", type=float, default=0)

aa = ap.parse_args()

fig = plt.figure(figsize=(4.2,3))

def get_phase(Ne,No):
	with open(f"output_Lz/{Ne}e{No}_Lz.dat") as f:
		data = [list(map(float, x.split())) for x in f.readlines()]
		theta2 = np.array([x[0] for x in data])
		Lz2   = np.array([-x[1] for x in data])

	Omega2  = (1-np.cos(theta2))/2 # This is actually Solid angle / 4π
	#Lz2_mod = (Lz2 - Lz2[0])*(No-1)/(No+1)

	m, b = np.polyfit(Omega2[-5:], Lz2[-5:],1) # fit y = mx + b for the last 4 points of data
	print(f"m = {m}")
	Lz2_mod = Lz2 - Lz2[0] - Omega2*(Ne/2) # Here Nϕ = Nₒ
	print(Lz2_mod)
	return Omega2, Lz2_mod

plt.xlabel(r"$\Omega/4\pi$")
plt.ylabel(r"Braiding phase/$2\pi$")


for Ne in [8,10,12]:
	Omega, phase = get_phase(Ne, 2*Ne-1)
	plt.plot(Omega, phase, label=f"{Ne} electrons")

if abs(aa.line)<2:
	ref = aa.line * np.ones(Omega.size)
	plt.plot(Omega, ref, "k--", label=f"{aa.line}")

plt.legend()


plt.savefig(f"plots/LZ_even.svg")
