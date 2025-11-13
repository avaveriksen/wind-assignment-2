import math
import cmath
import matplotlib.pyplot as plt

# --------------------------
# Machine / turbine params
# --------------------------
poles = 266
Ls = 7.8e-3
Rs = 0.254
psi_f = 19.5
P_mech_rated = 10e6
P_rot_loss_rated = 225e3
N_rated_rpm = 12.0


# --------------------------
# Shared turbine model
# --------------------------
def turbine_mechanics(Vw):
    N_rpm = Vw
    omega_m = 2*math.pi*N_rpm/60
    omega_rated = 2*math.pi*N_rated_rpm/60
    pole_pairs = poles/2
    omega_e = pole_pairs * omega_m

    Ea = omega_e * psi_f
    Xs = omega_e * Ls

    k_w = P_mech_rated / (omega_rated**3)
    P_mech_total = k_w*(omega_m**3)

    P_rot_loss = P_rot_loss_rated * (omega_m/omega_rated)**2
    P_mech_net = P_mech_total - P_rot_loss
    P_phase = P_mech_net / 3

    return dict(Ea=Ea, Xs=Xs, P_mech_net=P_mech_net, P_phase=P_phase)


# --------------------------
# PF = 1 mode
# --------------------------
def generator_pf1(Vw):
    m = turbine_mechanics(Vw)
    Ea=m["Ea"]; Xs=m["Xs"]; P_phase=m["P_phase"]; P_mech_net=m["P_mech_net"]

    a = Xs**2
    b = -Ea**2
    c = P_phase**2
    disc = b*b - 4*a*c
    y1 = (-b + math.sqrt(disc))/(2*a)
    y2 = (-b - math.sqrt(disc))/(2*a)

    candidates = []
    for y in (y1, y2):
        if y <= 0: continue
        Ia = math.sqrt(y)
        cosd = P_phase/(Ea*Ia)
        if abs(cosd) > 1: continue
        Va = (Ea*cosd - Ia*Rs) + 0j
        S = 3*Va*Ia
        candidates.append((Ia, Va, S, cosd))

    Ia, Va, S, cosd = max(candidates, key=lambda x:x[3])
    return dict(Vw=Vw, Ea=Ea, Va=Va, Ia=Ia, S=S, P_mech_net=P_mech_net)


# --------------------------
# β = 0 mode
# --------------------------
def generator_beta0(Vw):
    m = turbine_mechanics(Vw)
    Ea=m["Ea"]; Xs=m["Xs"]; P_mech_net=m["P_mech_net"]

    Ia = P_mech_net / (3*Ea)
    Va = Ea - Ia*(Rs + 1j*Xs)
    S = 3*Va*Ia

    return dict(Vw=Vw, Ea=Ea, Va=Va, Ia=Ia, S=S, P_mech_net=P_mech_net)


# -----------------------------------------------------
# PRINT TABLE 1: PF = 1
# -----------------------------------------------------
print("\n==========================")
print(" TABLE 1 — PF = 1 Mode")
print("==========================")
print(f"{'Vw':>4} | {'Ea[V]':>10} | {'Va_real[V]':>12} | {'Va_imag[V]':>12} | {'Ia[A]':>10} | {'S_real[kW]':>14} | {'S_imag[kVAr]':>14} | {'Pmech[kW]':>12}")
print("-"*110)

pf_list = []
for Vw in range(1,13):
    r = generator_pf1(Vw)
    pf_list.append(r)
    print(f"{Vw:4.0f} | {r['Ea']:10.2f} | {r['Va'].real:12.2f} | {r['Va'].imag:12.2f} | {r['Ia']:10.2f} | "
          f"{r['S'].real/1e3:14.2f} | {r['S'].imag/1e3:14.2f} | {r['P_mech_net']/1e3:12.2f}")


# -----------------------------------------------------
# PRINT TABLE 2: β = 0
# -----------------------------------------------------
print("\n==========================")
print(" TABLE 2 — Beta = 0 Mode")
print("==========================")
print(f"{'Vw':>4} | {'Ea[V]':>10} | {'Va_real[V]':>12} | {'Va_imag[V]':>12} | {'Ia[A]':>10} | {'S_real[kW]':>14} | {'S_imag[kVAr]':>14} | {'Pmech[kW]':>12}")
print("-"*110)

b0_list = []
for Vw in range(1,13):
    r = generator_beta0(Vw)
    b0_list.append(r)
    print(f"{Vw:4.0f} | {r['Ea']:10.2f} | {r['Va'].real:12.2f} | {r['Va'].imag:12.2f} | {r['Ia']:10.2f} | "
          f"{r['S'].real/1e3:14.2f} | {r['S'].imag/1e3:14.2f} | {r['P_mech_net']/1e3:12.2f}")


# -----------------------------------------------------
# PLOT Active Power vs Wind Speed
# -----------------------------------------------------
plt.figure(figsize=(8,5))

plt.plot([r["Vw"] for r in pf_list],
         [r["S"].real/1e6 for r in pf_list],
         marker="o", label="PF = 1")

plt.plot([r["Vw"] for r in b0_list],
         [r["S"].real/1e6 for r in b0_list],
         marker="s", label="Beta = 0")

plt.xlabel("Wind Speed (m/s)")
plt.ylabel("Active Power (MW)")
plt.title("Active Power vs Wind Speed for Both Modes")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
