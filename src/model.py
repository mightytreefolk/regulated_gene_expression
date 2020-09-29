
import math

n_A = 6.023E23 #Avogadro's Number

"""
This models the ODE of an unregulated genetic expression of a prokaryotic cell
k0 is the rate at which mRNA is being generated
k1 is the rate at which proteins are being 
dm and dp are the degradation constants of mRNA and proteins respectively
C1 and C2 are integration constants
"""


class UnregulatedGeneExpression:

    # t0 is a single number and t is a vector of numbers
    # const is an array of numbers that represent the constants of the reactions
    def __init__(self, t0, t, const):
        self.t0 = t0
        self.t = t
        self.k0 = const[0]
        self.k1 = const[1]
        self.dm = const[2]
        self.dp = const[3]

    def unregulated_m_rna(self):
        m_rna = C1 * math.exp(-self.dm * self.t0) + self.k0/self.dm
        return m_rna

    def unregulated_protein(self):
        protein = (self.k1 * self.k0) / (self.dm * self.dp) + (C1 * math.exp(-self.dm * self.t0))/(self.dp - self.dm) \
                  + C2 * math.exp(-self.dp * self.t0)
        return protein

    # Use scipy.odeint to evaluate
    def unregulated_gene_expression(self, z, t=.1):
        m0 = z[0]
        p0 = z[1]
        dmdt = self.k0 - self.dm * m0
        dpdt = self.k1 * m0 - self.dp*p0
        dzdt = dmdt, dpdt
        return dzdt

    # Use scipy.odeint to evaluate
    def reduced_model(self, z, t):
        p0 = z[1]
        dpdt = (self.k1 * self.k0) / self.dm - self.dp * p0
        return 0, dpdt

    def gillespie(self):
        a1 = self.k0
        a2 = self.k1 * Nm
        a3 = self.dm * Nm
        a4 = self.dp * Np

        total = a1 + a2 + a3 + a4







