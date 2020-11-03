import math
import random
import numpy
import pandas
import os
from datetime import datetime
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
        elif r[2] + r[1] + r[0] < n <= r[3] + r[2] + r[1] + r[0]:
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
    def __init__(self, tmax=10, m0=0, p0=0, const=None, number_of_sims=1, cell_cycle=6300):
        super().__init__(tmax, m0, p0, const)
        self.number_of_sims = number_of_sims
        self.two_genes = cell_cycle/2
        self.divide = cell_cycle

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
        clock = 0
        cell_cycle = [0]
        t = 0  # start time
        r0 = self.initial_state()
        g0 = 1
        genes.append(g0)
        mrna.append(r0[0])
        protein.append(r0[1])
        time.append(t)
        while t <= self.tmax:
            # Refactor this. Not doing what I want it to.
            if self.two_genes <= counter <= self.divide:
                genes.append(2)
                a = self.two_genes_propensities(r0[0], r0[1])
                next_rxn = self.next_reaction(a)
                sim_time = self.time_to_next_rxn(a, t)
                t = sim_time[0]
                if counter == self.divide:
                    r = self.division_of_molecules(r0)
                    mrna.append(r[0])
                    protein.append(r[1])
                    time.append(t)
                    r0 = r
                    # Mark relevant data
                    counter = 0
                    clock = 0
                    counter0.append(counter)
                    cell_cycle.append(clock)
                    divide.append("Yes")
                else:
                    r = self.update_reaction_vector(r0, next_rxn)
                    mrna.append(r[0])
                    protein.append(r[1])
                    time.append(t)
                    r0 = r
                    # Mark relevant data
                    counter = counter + 1
                    clock = sim_time[1] + clock
                    counter0.append(counter)
                    cell_cycle.append(clock)
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
                counter = counter + 1
                clock = sim_time[1] + clock
                counter0.append(counter)
                cell_cycle.append(clock)
                genes.append(1)
                divide.append("No")

        sim_df = pandas.DataFrame()
        sim_df["Seconds"] = time
        time = [x / 3600 for x in time]
        sim_df["Hours"] = time
        sim_df["Time"] = time
        sim_df["Proteins"] = protein
        sim_df["mRNA"] = mrna
        sim_df["Gene Number"] = genes
        sim_df["Divide"] = divide
        sim_df["Counter"] = counter0
        sim_df["Clock"] = cell_cycle
        return sim_df

    def multiple_cells(self):
        timestamp = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        directory = "Dir_with_runs-{time}".format(time=timestamp)
        os.mkdir(directory)
        for i in range(0, self.number_of_sims):
            path = os.path.join(directory, "Run{n}.csv".format(n=i))
            sim = self.sim()
            sim.to_csv(path, index=False, encoding='utf-8', sep='\t',)
        return directory


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

    def numerical_solution_1_gene(self, z, t):
        m0 = z[0]
        p0 = z[1]
        dmdt = self.k0 - self.dm * m0
        dpdt = self.k1 * m0 - self.dp * p0
        dzdt = dmdt, dpdt
        return dzdt

    def numerical_solution_2_gene(self, z, t):
        m0 = z[0]
        p0 = z[1]
        dmdt = 2 * self.k0 - self.dm * m0
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
        counter = 0
        counters = []
        counters.append(counter)
        # solve ODE
        for i in range(1, self.n):
            if counter < 1800:
                # span for next time step
                tspan = [self.t[i - 1], self.t[i]]

                # solve for next step
                z = odeint(self.numerical_solution_1_gene, z0, tspan)

                # store solution for plotting
                m[i] = z[1][0]
                p[i] = z[1][1]
                # next initial condition
                z0 = z[1]
                counter = counter + 1
                counters.append(counter)
            elif 1800 <= counter < 3600:
                tspan = [self.t[i - 1], self.t[i]]
                # solve for next step
                z = odeint(self.numerical_solution_2_gene, z0, tspan)
                # store solution for plotting
                m[i] = z[1][0]
                p[i] = z[1][1]
                # next initial condition
                z0 = z[1]
                counter = counter + 1
                counters.append(counter)
            elif counter >= 3600:
                m[i] = m[i-1]/2
                p[i] = p[i-1]/2
                z0 = m[i], p[i]
                counter = 0
                counters.append(counter)

        time = [x / 3600 for x in self.t]
        sim_df = pandas.DataFrame()
        sim_df["Counters"] = counters
        sim_df["Time"] = time
        sim_df["Proteins"] = p
        sim_df["mRNA"] = m
        return sim_df