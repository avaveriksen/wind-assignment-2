# Article on how to solve non-linear equations in Python: https://coderivers.org/blog/fsolve-python/

from scipy.optimize import fsolve
import numpy as np
import math, pprint

from Beta0 import Vw, turbine_mechanics

poles = 266
Ls = 7.8e-3
Rs = 0.254
psi_f = 19.5
P_mech_rated = 10e6
P_rot_loss_rated = 225e3
N_rated_rpm = 12.0
class Tools:
    def equations(self, VS, *args):
        '''
        Build equations for the solver of non-linear equations to use.
        :param VS: Input voltages V_S1 and V_S2
        :param args: Real Power input from Turbine 1 and 2
        :return: equations for V_B, f(V_S1, V_S2) - VB = 0 form
        '''
        V_T1 = VS[0] + 1j*VS[1]
        V_T2 = VS[2] + 1j*VS[3]

        P1, P2 = args

        V_POC = 4e3 / np.sqrt(3)

        Z11, Z12, Z21, Z22, Z31, Z32, Z41, Z42, Z5, Z6 = self.get_impedances(verb=False)

        # equations copied from /Anders/equations.ipynb
        eq1 = (-P1*Z11*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T2) + V_T1*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T1)*np.conj(V_T2) - Z21*Z32*(P1*Z11 + P1*Z31 - V_T1*np.conj(V_T1))*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T2) - Z41*Z42*(P1*Z11*Z32*Z5*np.conj(V_T2) + P1*Z31*Z32*Z5*np.conj(V_T2) + P2*Z12*Z31*Z5*np.conj(V_T1) + P2*Z31*Z32*Z5*np.conj(V_T1) + V_POC*Z31*Z32*np.conj(V_T1)*np.conj(V_T2) - V_T1*Z32*Z5*np.conj(V_T1)*np.conj(V_T2) - V_T2*Z31*Z5*np.conj(V_T1)*np.conj(V_T2)))/(Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T1)*np.conj(V_T2))
        eq2 = (-P2*Z12*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T1) + V_T2*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T1)*np.conj(V_T2) - Z22*Z31*(P2*Z12 + P2*Z32 - V_T2*np.conj(V_T2))*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T1) - Z41*Z42*(P1*Z11*Z32*Z5*np.conj(V_T2) + P1*Z31*Z32*Z5*np.conj(V_T2) + P2*Z12*Z31*Z5*np.conj(V_T1) + P2*Z31*Z32*Z5*np.conj(V_T1) + V_POC*Z31*Z32*np.conj(V_T1)*np.conj(V_T2) - V_T1*Z32*Z5*np.conj(V_T1)*np.conj(V_T2) - V_T2*Z31*Z5*np.conj(V_T1)*np.conj(V_T2)))/(Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5)*np.conj(V_T1)*np.conj(V_T2))

        return np.array([eq1.real, eq1.imag, eq2.real, eq2.imag])

    def get_impedances(self, verb):
        '''
        Calculate impedances of the system
        :param verbose: Boolean
        :return: Z
        '''

        w = 2 * np.pi * 50

        # 220kV export cable, referred twice
        r = 0.25
        l = 2.5e-3
        c = 0.05e-6
        len = 30

        Z_ex_shunt = 1 / (1j * w * c * len)
        Z_ex_ser = r * len + 1j * w * l * len

        # Wind farm transformer
        # secondary side has already been referred to primary, no need to do more
        l1 = 2154e-6 #Henry
        l2 = l1
        r1 = 834e-3 #Ohm
        r2 = r1
        Z_wft = r1 + r2 + 1j * w * (l1 + l2)

        # Turbine cables, 33kV
        r = 0.1 # Ohm pr.km
        l = 0.6e-3 # H pr.km
        c = 0.3e-6 # F pr.km
        len1 = 6
        len2 = 16

        Z_33k_shunt_1 = 1 / (1j * w * c * len1)
        Z_33k_shunt_2 = 1 / (1j * w * c * len2)
        Z_33k_ser_1 = len1 * (r + 1j * w * l)
        Z_33k_ser_2 = len2 * (r + 1j * w * l)

        # Turbine transformer
        # Secondary side has already been referred to primary, no need to do more
        l1 = 10e-6
        l2 = l1
        r1 = 42e-3
        r2 = r1
        # core and magnetization
        Rc = 55
        Lm = 55e-3

        Z_wtt_sec = r2 + 1j * w * l2
        Z_wtt_prim = r1 + 1j * w * l1
        Z_wtt_shunt = (Rc**-1 + (1j * w * Lm)**-1)**-1

        # Referring
        Z_ex_shunt = Z_ex_shunt * (33 / 220)**2 * (4 / 33)**2
        Z_ex_ser = Z_ex_ser * (33 / 220)**2 * (4 / 33)**2
        
        Z_wft = Z_wft * (4 / 33)**2
        
        Z_33k_ser_1 = Z_33k_ser_1 * (4 / 33)**2
        Z_33k_ser_2 = Z_33k_ser_2 * (4 / 33)**2
        Z_33k_shunt_1 = Z_33k_shunt_1 * (4 / 33)**2
        Z_33k_shunt_2 = Z_33k_shunt_2 * (4 / 33)**2

        # Equivalent impedances
        Z11 = Z_wtt_prim
        Z12 = Z_wtt_prim
        Z21 = Z_wtt_sec + Z_33k_ser_1
        Z22 = Z_wtt_sec + Z_33k_ser_2
        Z31 = Z_wtt_shunt
        Z32 = Z_wtt_shunt
        Z41 = Z_33k_shunt_1
        Z42 = Z_33k_shunt_2
        Z5 = Z_wft + Z_ex_ser
        Z6 = Z_ex_shunt

        Z = np.array([Z11, Z12, Z21, Z22, Z31, Z32, Z41, Z42, Z5, Z6])
        Z = np.round(Z, 5)
        if verb:
            print("---------- IMPEDANCES ----------")
            print("All referred to 4kV")

            print("\n----- Wind Turbine Transformer -----")
            print("Z_wtt_prim (Z1):", Z_wtt_prim)
            print("Z_wtt_shunt (Zm):", Z_wtt_shunt)
            print("Z_wtt_sec (Z2'):", Z_wtt_sec)

            print("\n----- 33kV Cables -----")
            print("Z_33k_ser_1 (Z_RL):", Z_33k_ser_1)
            print("Z_33k_shunt_1 (Z_C):", Z_33k_shunt_1)
            print("Z_33k_ser_2 (Z_RL):", Z_33k_ser_2)
            print("Z_33k_shunt_2 (Z_C):", Z_33k_shunt_2)

            print("\n----- Wind Farm Transformer -----")
            print("Z_wtt_prim (Zseries):", Z_wft)

            print("\n----- 220kV cable -----")
            print("Z_ex_ser (Z_RL):", Z_ex_ser)
            print("Z_ex_shunt (Z_C):", Z_ex_shunt)

            print("\n\n ---- LUMPED IMPEDANCES -----")

            Z_name = ['Z11', 'Z12', 'Z21', 'Z22', 'Z31', 'Z32', 'Z41', 'Z42', 'Z5', 'Z6']
            for index, impedance in enumerate(Z):
                print(Z_name[index] , ": " , impedance)
        return Z

    def wind_turbine_generator(self, Vw):
        """
        Calculate generator electrical output for a given wind speed (m/s).
        Returns P_total [MW], Va [V], Ia [A].
        """
        # --- Given parameters ---
        poles = 266
        Ls = 7.8e-3  # H
        Rs = 0.254  # ohm
        psi_f = 19.5  # Wb
        P_mech_rated = 10e6  # W
        P_rot_loss_rated = 225e3  # W at 12 rpm
        N_rated_rpm = 12.0  # rated mechanical speed

        # --- Mechanical speed relationship ---
        N_rpm = Vw  # assume 1:1 (RPM = windspeed)
        omega_m = 2 * math.pi * N_rpm / 60  # mechanical rad/s
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
    def generator_beta0(Vw):
        m = turbine_mechanics(Vw)
        Ea=m["Ea"]; Xs=m["Xs"]; P_mech_net=m["P_mech_net"]

        Ia = P_mech_net / (3*Ea)
        Va = Ea - Ia*(Rs + 1j*Xs)
        S = 3*Va*Ia

        return dict(Vw=Vw, Ea=Ea, Va=Va, Ia=Ia, S=S, P_mech_net=P_mech_net)