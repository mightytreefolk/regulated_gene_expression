import math
import random
import numpy
import pandas
from scipy.integrate import odeint


class Gillespie:

    # const is an array of numbers that represent the constants of the reactions
    def __init__(self, tmax=10, m0=0, p0=0, const=None):
        if const is None:
            const = [1, 1, .1, .1]
        self.tmax = tmax
        self.k0 = const[0]
        self.k1 = const[1]
        self.dm = const[2]
        self.dp = const[3]
        self.Nm = m0
        self.Np = p0

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
        return t, dt

    def update_reaction_vector(self, current_state, update_vector):
        current_state = numpy.array(current_state)
        update_vector = numpy.array(update_vector)
        updated = current_state + update_vector
        return updated


class CellDivision(Gillespie):
    def __init__(self, tmax=10, m0=0, p0=0, const=None):
        super().__init__(tmax, m0, p0, const)

    def two_genes_propensities(self, Nm, Np):
        a1 = 2 * self.k0
        a2 = self.k1 * Nm
        a3 = self.dm * Nm
        a4 = self.dp * Np

        total = a1 + a2 + a3 + a4

        r1 = a1 / total
        r2 = a2 / total
        r3 = a3 / total
        r4 = a4 / total

        return [r1, r2, r3, r4, total]

    def division_of_molecules(self, current_state):
        parent_cell = []
        current_mrna, current_protein = current_state[0], current_state[1]
        mrna = 0
        prot = 0
        for i in range(current_mrna):
            stays = bool(random.randint(0, 1))
            if stays == 0:
                mrna = mrna + 1
            else:
                pass
        parent_cell.append(mrna)
        for i in range(current_protein):
            stays = bool(random.randint(0, 1))
            if stays == 0:
                prot = prot + 1
            else:
                pass
        parent_cell.append(prot)
        return parent_cell

    def sim(self):
        time = []
        protein = []
        mrna = []
        genes = []
        divide = ["No"]
        counter0 = [0]
        counter = 0
        t = 0  # start time
        r0 = self.initial_state()
        g0 = 1
        genes.append(g0)
        mrna.append(r0[0])
        protein.append(r0[1])
        time.append(t)
        while t <= self.tmax:
            # Refactor this. Not doing what I want it to.
            if 1800 <= counter <= 3601:
                genes.append(2)
                a = self.two_genes_propensities(r0[0], r0[1])
                next_rxn = self.next_reaction(a)
                sim_time = self.time_to_next_rxn(a, t)
                t = sim_time[0]
                if 3599 <= counter <= 3601:
                    r = self.division_of_molecules(r0)
                    mrna.append(r[0])
                    protein.append(r[1])
                    time.append(t)
                    r0 = r
                    # Mark relevant data
                    counter = 0
                    counter0.append(counter)
                    divide.append("Yes")
                else:
                    r = self.update_reaction_vector(r0, next_rxn)
                    mrna.append(r[0])
                    protein.append(r[1])
                    time.append(t)
                    r0 = r
                    # Mark relevant data
                    counter = sim_time[1] + counter
                    counter0.append(counter)
                    divide.append("No")
            else:
                a = self.update_propensities(r0[0], r0[1])
                next_rxn = self.next_reaction(a)
                sim_time = self.time_to_next_rxn(a, t)
                t = sim_time[0]
                r = self.update_reaction_vector(r0, next_rxn)
                mrna.append(r[0])
                protein.append(r[1])
                time.append(t)
                r0 = r
                # Mark relevant data
                counter = sim_time[1] + counter
                counter0.append(counter)
                genes.append(1)
                divide.append("No")

        time = [x/3600 for x in time]

        sim_df = pandas.DataFrame()
        sim_df["Time"] = time
        sim_df["Proteins"] = protein
        sim_df["mRNA"] = mrna
        sim_df["Gene Number"] = genes
        sim_df["Divide"] = divide
        sim_df["Counter"] = counter0
        sim_df.to_csv("dataframe.csv", index=False, encoding='utf-8', sep='\t',)
        return sim_df


class DeterministicCellDivision:
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

    def analytical_solution(self, t, counter, current_mrna, current_prot):
        C1 = self.m0 - self.k0 / self.dm
        C2 = self.p0 - (self.k1 * self.k0) / (self.dm * self.dp) + C1 * self.k1 / (self.dm - self.dp)
        if math.floor(counter) < 1800:
            m_rna = C1 * math.exp(-self.dm * t) + self.k0 / self.dm
            protein = (self.k1 * self.k0) / (self.dm * self.dp) - (C1 * self.k1 * math.exp(-self.dm * t)) / (self.dm - self.dp) + C2 * math.exp(-self.dp * t)
            return m_rna, protein
        elif 1800 <= math.floor(counter) <= 3600:
            m_rna = C1 * math.exp(-self.dm * t) + 2 * self.k0 / self.dm
            protein = (2 * self.k1 * self.k0)/(self.dm * (self.dm - self.dp)) + \
                      (2 * self.k0 * self.k1)/(self.dp * (self.dm-self.dp)) - \
                      (C1 * self.k1 * math.exp(-self.dm * t))/(self.dm - self.dp) + \
                       C2 * math.exp(-self.dp * t)
            return m_rna, protein
        elif math.floor(counter) == 3601:
            return current_mrna/2, current_prot/2

    def analytical_simulation(self):
        # store solutions
        m = []
        p = []
        counters = []
        current_mrna = self.m0
        current_prot = self.p0

        counter = 0
        # iterate over time:
        for i in self.t:
            if counter <= 3600:
                counter = counter + 1
                z = self.analytical_solution(i, counter=counter, current_mrna=current_mrna, current_prot=current_prot)
                current_mrna = z[0]
                current_prot = z[1]
                m.append(z[0])
                p.append(z[1])
                counters.append(counter)

            elif counter > 3600:
                counter = 0
                z = self.analytical_solution(i, counter=counter, current_mrna=current_mrna, current_prot=current_prot)
                current_mrna = z[0]
                current_prot = z[1]
                m.append(z[0])
                p.append(z[1])
                counters.append(counter)
        time = [x/3600 for x in self.t]
        sim_df = pandas.DataFrame()
        sim_df["Time"] = time
        sim_df["Proteins"] = p
        sim_df["mRNA"] = m
        # sim_df["Counter"] = counters
        sim_df.to_csv("dataframe.csv", index=False, encoding='utf-8', sep='\t', )
        return sim_df

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

        sim_df = pandas.DataFrame()
        sim_df["Time"] = self.t
        sim_df["Proteins"] = p
        sim_df["mRNA"] = m
        return sim_df