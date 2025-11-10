import math

def wind_turbine_generator(Vw):
    """
    Calculate generator electrical output for a given wind speed (m/s).
    Returns P_total [MW], Va [V], Ia [A].
    """
    # --- Given parameters ---
    poles = 266
    Ls = 7.8e-3      # H
    Rs = 0.254       # ohm
    psi_f = 19.5     # Wb
    P_mech_rated = 10e6        # W
    P_rot_loss_rated = 225e3   # W at 12 rpm
    N_rated_rpm = 12.0         # rated mechanical speed

    # --- Mechanical speed relationship ---
    N_rpm = Vw                # assume 1:1 (RPM = windspeed)
    omega_m = 2 * math.pi * N_rpm / 60    # mechanical rad/s
    omega_rated = 2 * math.pi * N_rated_rpm / 60

    # --- Electrical angular speed ---
    pole_pairs = poles / 2
    omega_e = pole_pairs * omega_m

    # --- Electrical quantities ---
    Ea = omega_e * psi_f
    Xs = omega_e * Ls

    # --- Mechanical powers ---
    # total available power scales with cube of omega (P ∝ ω³)
    k_w = P_mech_rated / (omega_rated ** 3)
    P_mech_total = k_w * (omega_m ** 3)

    # rotating losses scale with square of speed (P_loss ∝ ω²)
    P_rot_loss = P_rot_loss_rated * (omega_m / omega_rated) ** 2
    P_mech_net = P_mech_total - P_rot_loss
    P_phase = P_mech_net / 3.0

    # --- Solve quadratic for Ia² ---
    a = Xs ** 2
    b = -Ea ** 2
    c = P_phase ** 2
    disc = b ** 2 - 4 * a * c
    if disc < 0:
        raise ValueError("No real solution for Ia².")
    y1 = (-b + math.sqrt(disc)) / (2 * a)
    y2 = (-b - math.sqrt(disc)) / (2 * a)

    candidates = []
    for y in (y1, y2):
        if y <= 0:
            continue
        Ia = math.sqrt(y)
        sin_delta = Ia * Xs / Ea
        if sin_delta > 1:
            continue
        cos_delta = math.sqrt(1 - sin_delta ** 2)
        delta_deg = math.degrees(math.asin(sin_delta))
        Va = Ea * cos_delta - Ia * Rs
        P_total = 3 * Va * Ia
        eff_by_net = P_total / P_mech_net * 100
        candidates.append({
            "Ia": Ia,
            "Va": Va,
            "P_total": P_total,
            "eff_by_net": eff_by_net,
            "delta_deg": delta_deg,
            "cos_delta": cos_delta
        })

    if not candidates:
        raise ValueError("No valid operating points found.")

    # Choose the physical solution (smaller delta, larger cosδ)
    chosen = max(candidates, key=lambda s: s["cos_delta"])

    # --- Return results ---
    return {
        "Vw": Vw,
        "Va": chosen["Va"],
        "Ia": chosen["Ia"],
        "P_total_MW": chosen["P_total"] / 1e6,
        "efficiency_%": chosen["eff_by_net"],
        "delta_deg": chosen["delta_deg"]
    }


# --- Test loop: from 1 to 12 m/s ---
print("=== Wind Turbine Generator Performance ===")
print(f"{'Vw (m/s)':>8} | {'P_total (MW)':>12} | {'Va (V)':>10} | {'Ia (A)':>10} | {'η (%)':>8} | {'δ (°)':>8}")
print("-" * 65)

results = []
for Vw in range(1, 13):
    res = wind_turbine_generator(Vw)
    results.append(res)
    print(f"{res['Vw']:8.1f} | {res['P_total_MW']:12.4f} | {res['Va']:10.2f} | {res['Ia']:10.2f} | "
          f"{res['efficiency_%']:8.2f} | {res['delta_deg']:8.2f}")

print("=============================================================")


# --- Plot: Power vs Wind Speed ---
import matplotlib.pyplot as plt

# Extract data from results
Vw_values = [r["Vw"] for r in results]
P_values = [r["P_total_MW"] for r in results]

# Create plot
plt.figure(figsize=(7, 5))
plt.plot(Vw_values, P_values, marker='o', linewidth=2)
plt.title("Electrical Power Output vs Wind Speed")
plt.xlabel("Wind Speed (m/s)")
plt.ylabel("Electrical Power Output (MW)")
plt.grid(True)
plt.tight_layout()
plt.show()
