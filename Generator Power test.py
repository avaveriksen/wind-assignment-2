import math, pprint

# Given
poles = 266
Ls = 7.8e-3      # H
Rs = 0.254       # ohm
psi_f = 19.5     # Wb
P_mech_rated = 10e6
P_rot_loss_rated = 225e3
N_rated_rpm = 12.0

# speeds
omega_m = 2*math.pi*N_rated_rpm/60
pole_pairs = poles/2
omega_e = pole_pairs * omega_m

# electrical
Ea = omega_e * psi_f
Xs = omega_e * Ls

# mechanical powers
P_mech_total = P_mech_rated
P_rot_loss = P_rot_loss_rated
P_mech_net = P_mech_total - P_rot_loss
P_phase = P_mech_net / 3.0

# quadratic in y = Ia^2 : a y^2 + b y + c = 0
a = Xs**2
b = -Ea**2
c = P_phase**2
disc = b**2 - 4*a*c
if disc < 0:
    raise ValueError("No real roots")

y1 = (-b + math.sqrt(disc)) / (2*a)
y2 = (-b - math.sqrt(disc)) / (2*a)

candidates = []
for y in (y1, y2):
    if y <= 0:
        continue
    Ia = math.sqrt(y)
    sin_delta = Ia * Xs / Ea
    if sin_delta > 1:           # discard unphysical
        continue
    cos_delta = math.sqrt(1 - sin_delta**2)
    delta_deg = math.degrees(math.asin(sin_delta))
    Va = Ea * cos_delta - Ia * Rs     # real terminal voltage per phase
    P_terminal_phase = Va * Ia
    P_total = 3 * P_terminal_phase
    eff_by_net = P_total / P_mech_net * 100
    candidates.append({
        "Ia": Ia,
        "sin_delta": sin_delta,
        "cos_delta": cos_delta,
        "delta_deg": delta_deg,
        "Va": Va,
        "P_total": P_total,
        "eff_by_net": eff_by_net
    })

# choose physical solution: larger cos_delta (smaller delta)
if not candidates:
    raise ValueError("No valid candidate solutions found.")
chosen = max(candidates, key=lambda s: s["cos_delta"])

# print chosen
print("=== Selected operating point ===")
print(f"Ea = {Ea:.2f} V, Xs = {Xs:.4f} Ω, Rs = {Rs:.3f} Ω")
print(f"Ia = {chosen['Ia']:.3f} A  ∠ 0°")
print(f"delta = {chosen['delta_deg']:.5f}°  (cosδ = {chosen['cos_delta']:.6f})")
print(f"Va = {chosen['Va']:.3f} V  ∠ 0°")
print(f"P_total (3-ph) = {chosen['P_total']/1e6:.6f} MW")
print(f"Efficiency (P_total / P_mech_net) = {chosen['eff_by_net']:.6f} %")
print("All candidate solutions (for inspection):")
pprint.pprint(candidates, width=120)

