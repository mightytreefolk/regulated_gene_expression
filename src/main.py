from model import UnregulatedGeneExpression, GillespieUnregulatedGeneExpression
import numpy as np
from scipy.integrate import odeint
import pandas
import plotly.express as px
import plotly.graph_objects as go
import random


n_A = 6.023E23  # Avogadro's Number
e_coli_vol = 6.5E-16  # Liters


def ode_sim():

    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 1, 0.0022, 0.00125]
    odemodel = UnregulatedGeneExpression(0, 0, const=const)

    # initial condition
    z0 = [0, 0]

    # number of time points
    n = 10000

    # time points
    t = np.linspace(0, 10000, n)

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

    dfp = pandas.DataFrame()
    dfp["Time"] = t
    dfp["Proteins"] = p

    dfm = pandas.DataFrame()
    dfm["Time"] = t
    dfm["mRNA"] = m
    return dfm, dfp


def Gillespie_sim():
    # k0 (mRNA), k1 (protein), dm, dp
    initial_conditions = [0.0167, 1, 0.0022, 0.00125]
    # m0 = 1.79e-9, p0 = 1.5e-5
    gill_model = GillespieUnregulatedGeneExpression(tmax=1000, m0=0, p0=0, const=initial_conditions)

    time = []
    protein = []
    mrna = []
    t = 0  # start time
    r0 = gill_model.initial_state()
    mrna.append(r0[0])
    protein.append(r0[1])
    time.append(t)

    while t <= 10000:

        a = gill_model.update_propensities(r0[0], r0[1])
        next_rxn = gill_model.next_reaction(a)
        t = gill_model.time_to_next_rxn(a, t)
        r = gill_model.update_reaction_vector(r0, next_rxn)
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



def main():
    # path = "/Users/pthompson/Dropbox/Documents/PhD/First_Year/Goldings_lab/Math-Model-Sys-Bio/plots/ode_plot.html"
    #
    # """Plot ODE simulation"""
    # mrna, proteins = ode_sim()
    # trace1 = go.Scatter(
    #     x=mrna["Time"],
    #     y=mrna["mRNA"],
    #     name="Number of mRNA",
    #     marker=dict(
    #         color='rgb(34,163,192)'
    #     )
    # )
    #
    # trace2 = go.Scatter(
    #     x=proteins["Time"],
    #     y=proteins['Proteins'],
    #     name="Number of Protein"
    # )
    #
    # fig1 = px.scatter()
    # fig1.add_trace(trace1)
    # fig1.add_trace(trace2)
    # fig1.update_layout(
    #     title="ODE comparison of mRNA and Protein molecules over time",
    #     xaxis_title="Time (s)",
    #     yaxis_title="Number of Molecules",
    #     legend_title="Legend",
    #     font=dict(
    #         family="Courier New, monospace",
    #         size=12,
    #         color="Black"
    #     )
    # )
    # fig1.show()
    # fig1.write_html(path, include_plotlyjs=True)



    """Extract Gillespie data"""
    gill_mrna, gill_protein = Gillespie_sim()
    gill_path = "/Users/pthompson/Dropbox/Documents/PhD/First_Year/Goldings_lab/Math-Model-Sys-Bio/plots/gill_plot2.html"
    ode_mrna, ode_proteins = ode_sim()
    """Plot Gillespie sim"""
    gill_trace1 = go.Scatter(
        x=gill_mrna["Time"],
        y=gill_mrna["mRNA"],
        name="GILL Number of mRNA",
        marker=dict(
            color='rgb(34,163,192)'
        )
    )

    gill_trace2 = go.Scatter(
        x=gill_protein["Time"],
        y=gill_protein['Proteins'],
        name="GILL Number of Protein"
    )

    ode_trace1 = go.Scatter(
        x=ode_mrna["Time"],
        y=ode_mrna["mRNA"],
        name="ODE Number of mRNA",
        marker=dict(
            color='rgb(34,163,192)'
        )
    )

    ode_trace2 = go.Scatter(
        x=ode_proteins["Time"],
        y=ode_proteins['Proteins'],
        name="ODE Number of Protein"
    )

    gill_fig = px.scatter()
    gill_fig.add_trace(gill_trace1)
    gill_fig.add_trace(gill_trace2)
    gill_fig.add_trace(ode_trace1)
    gill_fig.add_trace(ode_trace2)
    gill_fig.update_layout(
        title="Gillespie comparison of mRNA and Protein molecules over time",
        xaxis_title="Time (s)",
        yaxis_title="Number of Molecules",
        legend_title="Legend",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    gill_fig.show()
    gill_fig.write_html(gill_path)



if __name__ == '__main__':
    main()
