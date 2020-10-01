import math
import random
import numpy
"""
This models the ODE of an unregulated genetic expression of a prokaryotic cell
k0 is the rate at which mRNA is being generated
k1 is the rate at which proteins are being 
dm and dp are the degradation constants of mRNA and proteins respectively
C1 and C2 are integration constants
m0 and p0 are the mRNA and protein concentrations in M
"""


class UnregulatedGeneExpression:

    # t is a number
    # const is an array of numbers that represent the constants of the reactions
    def __init__(self, t=0, m0=0, p0=0, const=(1, 1, .1, .1)):
        self.t = t
        self.m0 = m0
        self.p0 = p0
        self.k0 = const[0]
        self.k1 = const[1]
        self.dm = const[2]
        self.dp = const[3]

    def solved_unregulated(self):
        C1 = self.m0 - self.k0 / self.dm
        C2 = self.p0 - (self.k1 * self.k0) / (self.dm * self.dp) - C1 / (self.dp - self.dm)
        m_rna = C1 * math.exp(-self.dm * self.t) + self.k0 / self.dm
        protein = (self.k1 * self.k0) / (self.dm * self.dp) + (C1 * math.exp(-self.dm * self.t)) / (self.dp - self.dm) \
                  + C2 * math.exp(-self.dp * self.t)
        return m_rna, protein

    # Use scipy.odeint to evaluate
    def unregulated_gene_expression(self, z, t):
        m0 = z[0]
        p0 = z[1]
        dmdt = self.k0 - self.dm * m0
        dpdt = self.k1 * m0 - self.dp * p0
        dzdt = dmdt, dpdt
        return dzdt

    # Use scipy.odeint to evaluate
    def reduced_model(self, z, t):
        p0 = z[1]
        dpdt = (self.k1 * self.k0) / self.dm - self.dp * p0
        return 0, dpdt


class GillespieUnregulatedGeneExpression:

    # const is an array of numbers that represent the constants of the reactions
    def __init__(self, tmax=10, m0=0, p0=0, const=None):
        if const is None:
            const = [1, 1, .1, .1]
        self.tmax = tmax
        self.k0 = const[0]
        self.k1 = const[1]
        self.dm = const[2]
        self.dp = const[3]
        self.Nm = m0 * 6.5E-16 * 6.023E23
        self.Np = p0 * 6.5E-16 * 6.023E23

    def initial_state(self):
        return [self.Nm, self.Np]

    def update_propensities(self, Nm, Np):
        a1 = self.k0
        a2 = self.k1 * Nm
        a3 = self.dm * Nm
        a4 = self.dp * Np

        total = a1 + a2 + a3 + a4

        r1 = a1/total
        r2 = a2/total
        r3 = a3/total
        r4 = a4/total

        return [r1, r2, r3, r4, total]

    def next_reaction(self, r):
        n = random.uniform(0, 1)
        if 0 <= n <= r[0]:
            return [1, 0]
        elif r[0] < n <= r[1] + r[0]:
            return [0, 1]
        elif r[1] + r[0] < n <= r[2] + r[1] + r[0]:
            return [-1, 0]
        elif r[2] + r[1] + r[0] < n <= 1:
            return [0, -1]

    def time_to_next_rxn(self, a, t):
        dt = -math.log(1 - random.uniform(0, 1))/a[4]
        if t == self.tmax:
            return
        elif t != self.tmax:
            t = t + dt
            return t

    def update_reaction_vector(self, current_state, update_vector):
        current_state = numpy.array(current_state)
        update_vector = numpy.array(update_vector)
        updated = current_state + update_vector
        return updated


