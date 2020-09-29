from model import UnregulatedGeneExpression
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random


# k0, k1, dm, dp
const = [0.1, 0.1, 0.1, 0.1]
model = UnregulatedGeneExpression(0, 0, const)


def main():
    # initial condition
    z0 = [0.1, 0.3]

    # number of time points
    n = 10000

    # time points
    t = np.linspace(0, 1000, n)

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
        z = odeint(model.unregulated_gene_expression, z0, tspan)

        # store solution for plotting
        m[i] = z[1][0]
        p[i] = z[1][1]

        # next initial condition
        z0 = z[1]
    plt.subplot(211)
    plt.plot(t, m, 'r--', label='A(t)')
    plt.grid(True)
    plt.legend(loc='best')

    plt.subplot(212)
    plt.plot(t, p, 'r--', label='R(t)')
    plt.grid(True)
    plt.legend(loc='best')

    plt.show()

if __name__ == '__main__':
    main()
