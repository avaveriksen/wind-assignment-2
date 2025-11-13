# Article on how to solve non-linear equations in Python: https://coderivers.org/blog/fsolve-python/

from scipy.optimize import fsolve
import numpy as np
import math, pprint


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

        V_POC = 4e3

        Z11, Z12, Z21, Z22, Z31, Z32, Z41, Z42, Z5, Z6 = self.get_impedances(verb=False)

        # equations copied from /Anders/equations.ipynb
        eq1 = (-P1*V_T2*Z11*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5) + V_T1**2*V_T2*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5) - V_T2*Z21*Z32*(P1*Z11 + P1*Z31 - V_T1**2)*(Z41*Z42 + Z41*Z5 + Z42*Z5) - Z41*Z42*(P1*V_T2*Z11*Z32*Z5 + P1*V_T2*Z31*Z32*Z5 + P2*V_T1*Z12*Z31*Z5 + P2*V_T1*Z31*Z32*Z5 + V_POC*V_T1*V_T2*Z31*Z32 - V_T1**2*V_T2*Z32*Z5 - V_T1*V_T2**2*Z31*Z5))/(V_T1*V_T2*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5))
        eq2 = (-P2*V_T1*Z12*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5) + V_T1*V_T2**2*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5) - V_T1*Z22*Z31*(P2*Z12 + P2*Z32 - V_T2**2)*(Z41*Z42 + Z41*Z5 + Z42*Z5) - Z41*Z42*(P1*V_T2*Z11*Z32*Z5 + P1*V_T2*Z31*Z32*Z5 + P2*V_T1*Z12*Z31*Z5 + P2*V_T1*Z31*Z32*Z5 + V_POC*V_T1*V_T2*Z31*Z32 - V_T1**2*V_T2*Z32*Z5 - V_T1*V_T2**2*Z31*Z5))/(V_T1*V_T2*Z31*Z32*(Z41*Z42 + Z41*Z5 + Z42*Z5))

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
        Z_ex_ser = r * l + 1j * w * l * len

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

        Z_30k_shunt_1 = 1 / (1j * w * c * len1)
        Z_30k_shunt_2 = 1 / (1j * w * c * len2)
        Z_30k_ser_1 = len1 * (r + 1j * w * l)
        Z_30k_ser_2 = len2 * (r + 1j * w * l)

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
        Z_30k_ser_1 = Z_30k_ser_1 * (4 / 33)**2
        Z_30k_ser_2 = Z_30k_ser_2 * (4 / 33)**2
        Z_30k_shunt_1 = Z_30k_shunt_1 * (4 / 33)**2
        Z_30k_shunt_2 = Z_30k_shunt_2 * (4 / 33)**2

        # Equivalent impedances
        Z11 = Z_wtt_prim
        Z12 = Z_wtt_prim
        Z21 = Z_wtt_sec + Z_30k_ser_1
        Z22 = Z_wtt_sec + Z_30k_ser_2
        Z31 = Z_wtt_shunt
        Z32 = Z_wtt_shunt
        Z41 = Z_30k_shunt_1
        Z42 = Z_30k_shunt_2
        Z5 = Z_wft + Z_ex_ser
        Z6 = Z_ex_shunt

        Z = np.array([Z11, Z12, Z21, Z22, Z31, Z32, Z41, Z42, Z5, Z6])
        Z = np.round(Z, 5)
        if verb:
            Z_name = ['Z11', 'Z12', 'Z21', 'Z22', 'Z31', 'Z32', 'Z41', 'Z42', 'Z5', 'Z6']
            for index, impedance in enumerate(Z):
                print(Z_name[index] , ": " , impedance)
        return Z
