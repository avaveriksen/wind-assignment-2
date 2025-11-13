# All_impedances.py
# Compute impedances discussed (50 Hz). Prints rectangular (a + j b) and polar (|Z| ∠ θ).
# Copy/paste and run with Python 3.

import math, cmath

# ----- helpers -----
def to_polar(z):
    mag = abs(z)
    ang = math.degrees(cmath.phase(z))
    return mag, ang

def fmt_rect(z):
    return f"{z.real:.6g} {'+' if z.imag>=0 else '-'} j{abs(z.imag):.6g}"

def fmt_polar(z):
    mag, ang = to_polar(z)
    return f"{mag:.6g} ∠ {ang:.4g}°"

# ----- common -----
f = 50.0
omega = 2 * math.pi * f
print(f"Frequency: f = {f} Hz, omega = {omega:.6g} rad/s\n")

# ------------------------------------------------------------------
# 1) Wind turbine transformer (per-phase) - already referred to 4 kV
# Given (from earlier): R1 = R2' = 0.042 Ω ; L1 = L2' = 10e-6 H ; Lm = 55e-3 H ; Rc = 55 Ω
R1 = 0.042
R2p = 0.042
L1 = 10e-6
L2p = 10e-6
Lm = 55e-3
Rc = 55.0

Z1 = R1 + 1j * omega * L1
Z2p = R2p + 1j * omega * L2p
# magnetizing branch as Rc || jωLm
Zm = 1.0 / (1.0/Rc + 1.0/(1j * omega * Lm))

print("1) Wind turbine transformer (referred to 4 kV):")
print("   Given: R1=R2'=", R1, "Ω ; L1=L2'=", L1, "H ; Lm=", Lm, "H ; Rc=", Rc, "Ω")
print("   Z1 =", fmt_rect(Z1), " ->", fmt_polar(Z1))
print("   Z2' =", fmt_rect(Z2p), " ->", fmt_polar(Z2p))
print("   Zm (Rc || jωLm) =", fmt_rect(Zm), " ->", fmt_polar(Zm))
print("   Zseries = Z1 + Z2' =", fmt_rect(Z1+Z2p), " ->", fmt_polar(Z1+Z2p))
print()

# ------------------------------------------------------------------
# 2) 33 kV collection cables -> refer to 4 kV
# per-km: r = 0.1 Ω/km, l = 0.6 mH/km, c = 0.3 μF/km
r_col = 0.1          # Ω/km
l_col = 0.6e-3       # H/km
c_col = 0.3e-6       # F/km

a_4_33 = 4.0 / 33.0
a2_4_33 = a_4_33 ** 2

def cable_impedance_per_length(r_per_km, l_per_km, c_per_km, length_km):
    R = r_per_km * length_km
    L = l_per_km * length_km
    C = c_per_km * length_km
    Z_RL = R + 1j * omega * L
    Z_C = 1.0 / (1j * omega * C)
    return R, L, C, Z_RL, Z_C

collection_lengths = {"Turbine1": 6.0, "Turbine2": 16.0}

print("2) Collection cables (33 kV parameters) and referred to 4 kV:")
for name, Lkm in collection_lengths.items():
    R, Ltot, Ctot, Z_RL_33, Z_C_33 = cable_impedance_per_length(r_col, l_col, c_col, Lkm)
    Z_RL_4 = a2_4_33 * Z_RL_33
    Z_C_4  = a2_4_33 * Z_C_33
    print(f"  {name} (length {Lkm} km):")
    print(f"    totals (33 kV side): R={R:.6g} Ω, L={Ltot:.6g} H, C={Ctot:.6g} F")
    print(f"    Z_RL (33 kV) = {fmt_rect(Z_RL_33)} -> {fmt_polar(Z_RL_33)}")
    print(f"    Z_C  (33 kV) = {fmt_rect(Z_C_33)} -> {fmt_polar(Z_C_33)}")
    print(f"    refer to 4 kV: multiply by a^2 where a=4/33, a^2={a2_4_33:.8g}")
    print(f"    Z_RL (4 kV) = {fmt_rect(Z_RL_4)} -> {fmt_polar(Z_RL_4)}")
    print(f"    Z_C  (4 kV) = {fmt_rect(Z_C_4)} -> {fmt_polar(Z_C_4)}")
    print()

# ------------------------------------------------------------------
# 3) Wind farm transformer (33 kV / 220 kV) given referred to 33 kV
# Given: R1=R2'=0.834 Ω ; L1=L2'=2154e-6 H (neglect shunt)
R1_wf = 0.834
R2p_wf = 0.834
L1_wf = 2154e-6
L2p_wf = 2154e-6

Z1_33 = R1_wf + 1j * omega * L1_wf
Z2p_33 = R2p_wf + 1j * omega * L2p_wf
Zseries_33 = Z1_33 + Z2p_33

Z1_4 = a2_4_33 * Z1_33
Z2p_4 = a2_4_33 * Z2p_33
Zseries_4 = a2_4_33 * Zseries_33

print("3) Wind farm transformer (33kV side values) and referred to 4 kV:")
print("   Given (33 kV side): R1=R2'=", R1_wf, "Ω ; L1=L2'=", L1_wf, "H")
print("   Z1 (33kV) =", fmt_rect(Z1_33), "->", fmt_polar(Z1_33))
print("   Z2' (33kV) =", fmt_rect(Z2p_33), "->", fmt_polar(Z2p_33))
print("   Zseries (33kV) =", fmt_rect(Zseries_33), "->", fmt_polar(Zseries_33))
print("   Refer to 4kV by multiplying by a^2 (a=4/33):")
print("   Z1 (4kV) =", fmt_rect(Z1_4), "->", fmt_polar(Z1_4))
print("   Z2' (4kV) =", fmt_rect(Z2p_4), "->", fmt_polar(Z2p_4))
print("   Zseries (4kV) =", fmt_rect(Zseries_4), "->", fmt_polar(Zseries_4))
print()

# ------------------------------------------------------------------
# 4) Export cable (on 220 kV side) -> refer to 33 kV and 4 kV
r_exp = 0.25        # Ω/km
l_exp = 2.5e-3      # H/km
c_exp = 0.05e-6     # F/km
L_exp_km = 30.0

R_exp = r_exp * L_exp_km
L_exp = l_exp * L_exp_km
C_exp = c_exp * L_exp_km

Z_RL_220 = R_exp + 1j * omega * L_exp
Z_C_220 = 1.0 / (1j * omega * C_exp)

# 220 -> 33 factor
a_33_220 = 33.0 / 220.0
a2_33_220 = a_33_220 ** 2

Z_RL_33_from220 = a2_33_220 * Z_RL_220
Z_C_33_from220 = a2_33_220 * Z_C_220

# 33 -> 4 factor (we already have)
Z_RL_4_from220 = a2_4_33 * Z_RL_33_from220
Z_C_4_from220 = a2_4_33 * Z_C_33_from220

print("4) Export cable (30 km) on 220 kV side and referred to 33 kV and 4 kV:")
print(f"   per-km: r={r_exp} Ω/km, l={l_exp} H/km, c={c_exp} F/km")
print(f"   totals (220kV): R={R_exp:.6g} Ω, L={L_exp:.6g} H, C={C_exp:.6g} F")
print("   Z_RL (220kV) =", fmt_rect(Z_RL_220), "->", fmt_polar(Z_RL_220))
print("   Z_C  (220kV) =", fmt_rect(Z_C_220), "->", fmt_polar(Z_C_220))
print("   220->33 scale factor (33/220)^2 =", a2_33_220)
print("   Z_RL (33kV) =", fmt_rect(Z_RL_33_from220), "->", fmt_polar(Z_RL_33_from220))
print("   Z_C  (33kV) =", fmt_rect(Z_C_33_from220), "->", fmt_polar(Z_C_33_from220))
print("   33->4 scale factor (4/33)^2 =", a2_4_33)
print("   Z_RL (4kV) =", fmt_rect(Z_RL_4_from220), "->", fmt_polar(Z_RL_4_from220))
print("   Z_C  (4kV) =", fmt_rect(Z_C_4_from220), "->", fmt_polar(Z_C_4_from220))
print()

# ------------------------------------------------------------------
# If you want the results programmatically, you can assign them to a dict (uncomment if needed)
# results = {
#     "wind_turbine_transformer": {"Z1": Z1, "Z2p": Z2p, "Zm": Zm, "Zseries": Z1+Z2p},
#     "collection_cables": { name: {"Z_RL_33": cable_impedance_per_length(r_col,l_col,c_col,Lkm)[3],
#                                   "Z_C_33": cable_impedance_per_length(r_col,l_col,c_col,Lkm)[4],
#                                   "Z_RL_4": a2_4_33*cable_impedance_per_length(r_col,l_col,c_col,Lkm)[3],
#                                   "Z_C_4": a2_4_33*cable_impedance_per_length(r_col,l_col,c_col,Lkm)[4]} for name,Lkm in collection_lengths.items()},
#     "wind_farm_transformer": {"Z1_33": Z1_33, "Z2p_33": Z2p_33, "Zseries_33": Zseries_33, "Z1_4": Z1_4, "Z2p_4": Z2p_4, "Zseries_4": Zseries_4},
#     "export_cable": {"Z_RL_220": Z_RL_220, "Z_C_220": Z_C_220, "Z_RL_33": Z_RL_33_from220, "Z_C_33": Z_C_33_from220, "Z_RL_4": Z_RL_4_from220, "Z_C_4": Z_C_4_from220}
# }
#
# print("If you need 'results', uncomment the dict above and use it in your script.")
