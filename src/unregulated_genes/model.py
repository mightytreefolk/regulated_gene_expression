import math
import random
import numpy
import pandas
from scipy.integrate import odeint

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
    def __init__(self, tmax=10, num_of_datapoints=1000, m0=0, p0=0, const=(1, 1, .1, .1)):
        self.n = num_of_datapoints
        self.tmax = tmax
        self.m0 = m0
        self.p0 = p0
        self.k0 = const[0]
        self.k1 = const[1]
        self.dm = const[2]
        self.dp = const[3]
        self.t = numpy.linspace(0, self.tmax, self.n)  # time points

    def analytical_solution(self, t):
        C1 = self.m0 - self.k0 / self.dm
        C2 = self.p0 - (self.k1 * self.k0) / (self.dm * self.dp) + C1 * self.k1 / (self.dm - self.dp)
        m_rna = C1 * math.exp(-self.dm * t) + self.k0 / self.dm
        protein = (self.k1 * self.k0) / (self.dm * self.dp) - (C1 * self.k1 * math.exp(-self.dm * t)) / (self.dm - self.dp) + C2 * math.exp(-self.dp * t)
        return m_rna, protein

    def analytical_sim(self):
        # store solutions
        m = numpy.empty_like(self.t)
        p = numpy.empty_like(self.t)

        # record initial conditions
        m[0] = self.m0
        p[0] = self.p0
        # iterate over time:
        for i in range(1, self.n):
            z = self.analytical_solution(i)

            # store solution for plotting
            m[i] = z[0]
            p[i] = z[1]

        dfp = pandas.DataFrame()
        dfp["Time"] = self.t
        dfp["Proteins"] = p

        dfm = pandas.DataFrame()
        dfm["Time"] = self.t
        dfm["mRNA"] = m
        return dfm, dfp

    # Use scipy.odeint to evaluate
    def numerical_solution(self, z, t):
        m0 = z[0]
        p0 = z[1]
        dmdt = self.k0 - self.dm * m0
        dpdt = self.k1 * m0 - self.dp * p0
        dzdt = dmdt, dpdt
        return dzdt

    def numerical_sim(self):
        # store solutions
        m = numpy.empty_like(self.t)
        p = numpy.empty_like(self.t)
        # record initial conditions
        z0 = [self.m0, self.p0]
        m[0] = self.m0
        p[0] = self.p0
        # solve ODE
        for i in range(1, self.n):
            # span for next time step
            tspan = [self.t[i - 1], self.t[i]]

            # solve for next step
            z = odeint(self.numerical_solution, z0, tspan)

            # store solution for plotting
            m[i] = z[1][0]
            p[i] = z[1][1]
            # next initial condition
            z0 = z[1]

        dfp = pandas.DataFrame()
        dfp["Time"] = self.t
        dfp["Proteins"] = p

        dfm = pandas.DataFrame()
        dfm["Time"] = self.t
        dfm["mRNA"] = m
        return dfm, dfp



class GillespieUnregulatedGeneExpression:

    # const is an array of numbers that represent the constants of the reactions
    def __init__(self, tmax=10, m0=0, p0=0, const=None, num_cells=1):
        if const is None:
            const = [1, 1, .1, .1]
        self.tmax = tmax
        self.k0 = const[0]
        self.k1 = const[1]
        self.dm = const[2]
        self.dp = const[3]
        self.Nm = m0
        self.Np = p0
        self.num_cells = num_cells

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
        # Create an mRNA
        if 0 <= n <= r[0]:
            return [1, 0]
        # Create a protein
        elif r[0] < n <= r[1] + r[0]:
            return [0, 1]
        # Delete an mRNA
        elif r[1] + r[0] < n <= r[2] + r[1] + r[0]:
            return [-1, 0]
        # Delete a protein
        elif r[2] + r[1] + r[0] < n <= 1:
            return [0, -1]
        # no reaction occurs
        elif r[3] + r[2] + r[1] + r[0] < n < 1:
            return [0, 0]

    def time_to_next_rxn(self, a, t):
        dt = -math.log(1 - random.uniform(0, 1))/a[4]
        t = t + dt
        return t

    def update_reaction_vector(self, current_state, update_vector):
        current_state = numpy.array(current_state)
        update_vector = numpy.array(update_vector)
        updated = current_state + update_vector
        return updated

    def run_sim(self):
        time = []
        protein = []
        mrna = []
        t = 0  # start time
        r0 = self.initial_state()
        mrna.append(r0[0])
        protein.append(r0[1])
        time.append(t)
        while t <= self.tmax:
            a = self.update_propensities(r0[0], r0[1])
            next_rxn = self.next_reaction(a)
            t = self.time_to_next_rxn(a, t)
            r = self.update_reaction_vector(r0, next_rxn)
            mrna.append(r[0])
            protein.append(r[1])
            time.append(t)
            r0 = r

        dfp = pandas.DataFrame()
        dfp["Time"] = time
        dfp["Proteins"] = protein

        dfm = pandas.DataFrame()
        dfm["Time"] = time
        dfm["mRNA"] = mrna
        return dfm, dfp

    def multiple_cells_sim(self):
        dfp_multiple = pandas.DataFrame()
        dfm_multiple = pandas.DataFrame()
        dfmt_multiple = pandas.DataFrame()
        dfpt_multiple = pandas.DataFrame()
        if self.num_cells == 1:
            return self.run_sim()
        elif self.num_cells > 1:
            for i in range(0, self.num_cells):
                mrna, prot = self.run_sim()
                dfm_multiple["Run{n}".format(n=i)] = mrna["mRNA"]
                dfp_multiple["Run{n}".format(n=i)] = prot["Proteins"]
                dfmt_multiple["mRNA_Run_time{t}".format(t=i)] = mrna["Time"]
                dfpt_multiple["prot_Run_time{t}".format(t=i)] = prot["Time"]

                dfm_multiple = dfm_multiple.dropna()
                dfp_multiple = dfp_multiple.dropna()
                dfmt_multiple = dfmt_multiple.dropna()
                dfpt_multiple = dfpt_multiple.dropna()


            # Add average col to dataframes
            dfm_multiple["Average"], dfmt_multiple["Average_Time"] = dfm_multiple.mean(axis=1), \
                                                                     dfmt_multiple.mean(axis=1)
            dfp_multiple["Average"], dfpt_multiple["Average_Time"] = dfp_multiple.mean(axis=1), \
                                                                     dfpt_multiple.mean(axis=1)

            # Join the DataFrams
            dfm_final = pandas.concat([dfm_multiple, dfmt_multiple], axis=1, sort=False)
            dfp_final = pandas.concat([dfp_multiple, dfpt_multiple], axis=1, sort=False)
            return dfm_final, dfp_final


