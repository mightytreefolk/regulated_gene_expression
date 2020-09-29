from model import UnregulatedGeneExpression, GillespieUnregulatedGeneExpression
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas
import random


n_A = 6.023E23  # Avogadro's Number

# k0 (mRNA), k1 (protein), dm, dp
const = [1, 1, .2, .01]
odemodel = UnregulatedGeneExpression(0, 0, const=const)

e_coli_vol = 6.5E-16  # Liters


def ode_sim():
    # initial condition
    z0 = [6e-6, 5e-7]

    # number of time points
    n = 10000

    # time points
    t = np.linspace(0, 10, n)

    # store solutions
    m = np.empty_like(t)
    p = np.empty_like(t)

    # record initial conditions
    m[0] = z0[0]
    p[0] = z0[1]

    # solve ODE
    for i in range(1, n):
        # span for next time step
        tspan = [t[i - 1], t[i]]

        # solve for next step
        z = odeint(odemodel.unregulated_gene_expression, z0, tspan)

        # store solution for plotting
        m[i] = z[1][0]
        p[i] = z[1][1]

        # next initial condition
        z0 = z[1]
    return m, p, t


def Gillespie_sim():
    # k0 (mRNA), k1 (protein), dm, dp
    initial_conditions = [1, 1, .2, .01]
    gill_model = GillespieUnregulatedGeneExpression(tmax=10, m0=6e-6, p0=5e-7, const=initial_conditions)

    data = []
    t = 0  # start time
    r0 = gill_model.initial_state()
    data.append([r0, t])

    while t <= 10:

        a = gill_model.update_propensities(r0[0], r0[1])
        next_rxn = gill_model.next_reaction(a)
        t = gill_model.time_to_next_rxn(a, t)
        r = gill_model.update_reaction_vector(r0, next_rxn)
        data.append([r, t])
        r0 = r

    df = pandas.DataFrame(data)
    return df



def main():

    """Plot ODE simulation"""
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.plot(ode_sim()[2], ode_sim()[0], 'r--', label='m(t)')
    ax1.plot(ode_sim()[2], ode_sim()[1], 'b--', label='p(t)')

    ax1.grid(True)
    ax1.legend(loc='best')
    ax1.set_title('Concentrations of mRNA and protein over time (ODE)')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Concentration')


    """Extract Gillespie data"""
    gill_data = Gillespie_sim()
    proteins = []
    mrna = []
    molecular_data = gill_data[0]
    for i in molecular_data:
        mrna.append(i[0])
        proteins.append(i[1])
    time = gill_data[1]

    """Plot Gillespie sim"""
    fig1 = plt.figure()
    ax2 = fig1.add_subplot()
    ax2.plot(time, mrna, 'r', label='m(t)')
    ax2.plot(time, proteins, 'b', label='p(t)')
    ax2.grid(True)
    ax2.legend(loc='best')
    ax2.set_title('Numhber of mRNA and protein over time (Gillespie)')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Number of molecules')
    plt.show()




if __name__ == '__main__':
    main()
