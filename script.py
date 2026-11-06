import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters from Table 1
k1 = 1.0
km1 = 0.01
k2 = 0.5
km2 = 0.2

kcat = 0.15
KM = 1.0

k3 = 0.01
k4 = 0.1
k5 = 0.008

# Initial concentrations
L0 = 0.3
R0 = 1.0
TFc0 = 1.0

# Species order:
# [L, R, LR, LRa, TFc, pTFc, pTFn, TFn]
y0 = [L0, R0, 0, 0, TFc0, 0, 0, 0]

# ODE system
def plant_growth(t, y):

    L, R, LR, LRa, TFc, pTFc, pTFn, TFn = y

    # Ligand and receptor
    dL = -k1*L*R + km1*LR + km2*LRa

    dR = -k1*L*R + km1*LR + km2*LRa

    dLR = k1*L*R - km1*LR - k2*LR

    dLRa = k2*LR - km2*LRa

    # Michaelis-Menten phosphorylation
    v = kcat * LRa * TFc / (KM + TFc)

    # Transcription factor dynamics
    dTFc = -v + k5*TFn

    dpTFc = v - k3*pTFc

    dpTFn = k3*pTFc - k4*pTFn

    dTFn = k4*pTFn - k5*TFn

    return [dL, dR, dLR, dLRa,
            dTFc, dpTFc, dpTFn, dTFn]


# Solve the system
t_eval = np.linspace(0, 500, 1000)

sol = solve_ivp(
    plant_growth,
    [0, 500],
    y0,
    t_eval=t_eval,
    rtol=1e-9,
    atol=1e-11
)

# Plot each variable
species = ["L", "R", "LR", "LRa",
           "TFc", "pTFc", "pTFn", "TFn"]

fig, axes = plt.subplots(4, 2, figsize=(11, 11))

for i, ax in enumerate(axes.flat):
    ax.plot(sol.t, sol.y[i])
    ax.set_title(species[i])
    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration")
    ax.grid(True)

plt.tight_layout()
plt.show()


# Mass conservation
L, R, LR, LRa, TFc, pTFc, pTFn, TFn = sol.y

L_total = L + LR + LRa
R_total = R + LR + LRa
TF_total = TFc + pTFc + pTFn + TFn

plt.figure(figsize=(9, 5))

plt.plot(sol.t, L_total, label="Total ligand")
plt.plot(sol.t, R_total, label="Total receptor")
plt.plot(sol.t, TF_total, label="Total TF")

plt.xlabel("Time")
plt.ylabel("Total concentration")
plt.title("Mass conservation")
plt.legend()
plt.grid(True)
plt.show()