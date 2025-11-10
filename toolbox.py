# Article on how to solve non-linear equations in Python: https://coderivers.org/blog/fsolve-python/

from scipy.optimize import fsolve
import numpy as np
import math, pprint


class Tools:
    def equations(VS, *P):
        '''
        Build equations for the solver of non-linear equations to use.
        :param vars: Input voltages V_S1 and V_S2
        :param P1: Real Power input from Turbine 1
        :param P2: Real Power input from Turbine 2
        :return: equations for V_B, f(V_S1, V_S2) - VB = 0 form
        '''

        V_S1 = VS[0]
        V_S2 = VS[1]

        P1 = P[0]
        P2 = P[1]

        Z11, Z12, Z21, Z22, Z31, Z32, Z41, Z42, Z5, Z6 = self.get_impedances()

        # equations copied from /Anders/equations.ipynb
        eq1 = (-P1*Z11*Z31 + V_S1**2*Z31 - Z21*(P1*Z11 + P1*Z31 - V_S1**2))/(V_S1*Z31) - V_B
        eq2 = (-P2*Z12*Z32 + V_S2**2*Z32 - Z22*(P2*Z12 + P2*Z32 - V_S2**2))/(V_S2*Z32)
        eq3 = Z41*Z42*(P1*V_S2*Z11*Z22*Z5 + P1*V_S2*Z22*Z31*Z5 + P2*V_S1*Z12*Z31*Z5 + P2*V_S1*Z22*Z31*Z5 + V_POC*V_S1*V_S2*Z22*Z31 - V_S1**2*V_S2*Z22*Z5 - V_S1*V_S2**2*Z31*Z5)/(V_S1*V_S2*Z22*Z31*(Z41*Z42 + Z41*Z5 + Z42*Z5))

        return [eq1, eq2, eq3]

    def get_impedances(self):
        '''
        Calculate impedances of the system
        :return: Z
        '''

        # Calculate all impedances
        Z11 = 1
        Z12 = 1
        Z21 = 1
        Z22 = 1
        Z31 = 1
        Z32 = 1
        Z41 = 1
        Z42 = 1
        Z5 = 1
        Z6 = 1

        return Z11, Z12, Z21, Z22, Z31, Z32, Z41, Z42, Z5, Z6
