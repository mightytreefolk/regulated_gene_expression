import numpy as np
import os
import pandas
import plotly.express as px
import plotly.graph_objects as go

from plotly.subplots import make_subplots
from datetime import datetime
from scipy.integrate import odeint
from model import UnregulatedGeneExpression, GillespieUnregulatedGeneExpression

n_A = 6.023E23  # Avogadro's Number
e_coli_vol = 6.5E-16  # Liters


def ode_sim(tmax, z0, n, const):
    odemodel = UnregulatedGeneExpression(0, 0, const=const)

    # time points
    t = np.linspace(0, tmax, n)

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


def numerical_sim(tmax, z0, n, const):
    numerical_model = UnregulatedGeneExpression(m0=0, p0=0, const=const)

    # time points
    t = np.linspace(0, tmax, n)

    # store solutions
    m = np.empty_like(t)
    p = np.empty_like(t)

    # record initial conditions
    m[0] = z0[0]
    p[0] = z0[1]

    # iterate over time:
    for i in range(1, n):
        z = numerical_model.solved_unregulated(i)

        # store solution for plotting
        m[i] = z[0]
        p[i] = z[1]

    dfp = pandas.DataFrame()
    dfp["Time"] = t
    dfp["Proteins"] = p

    dfm = pandas.DataFrame()
    dfm["Time"] = t
    dfm["mRNA"] = m
    return dfm, dfp


def main():
    # seconds for sims
    tmax = 10000
    # number of data points
    n = 10000
    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 1, 0.0022, 0.00125]
    # m0, p0
    initial_conditions = [7, 6072]

    """Extract data"""
    num_mrna, num_protein = numerical_sim(tmax=tmax, z0=initial_conditions, n=n, const=const)
    ode_mrna, ode_proteins = ode_sim(tmax=tmax, z0=initial_conditions, n=n, const=const)
    gill_model = GillespieUnregulatedGeneExpression(tmax=tmax, m0=initial_conditions[0], p0=initial_conditions[1],
                                                    const=const)
    gill_mrna, gill_protein = gill_model.run_sim()

    """Define plot paths"""
    base_path = "/Users/pthompson/Dropbox/Documents/PhD/First_Year/Goldings_lab/Math-Model-Sys-Bio/plots/"
    run_name = "steady_state"
    timestamp = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    num_ode_path = os.path.join(base_path, "ode_num_compare_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                                run_name=run_name))
    gill_path = os.path.join(base_path, "gill_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                  run_name=run_name))
    stat_path = os.path.join(base_path, "stat_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                  run_name=run_name))
    num_ode_image_path = os.path.join(base_path,
                                      "images/ode_num_compare_plot-{run_name}-{time}.png".format(time=timestamp,
                                                                                                 run_name=run_name))
    gill_image_path = os.path.join(base_path, "images/gill_plot-{run_name}-{time}.png".format(time=timestamp,
                                                                                              run_name=run_name))
    stat_image_path = os.path.join(base_path, "images/stat_plot-{run_name}-{time}.png".format(time=timestamp,
                                                                                              run_name=run_name))

    """Create Gillespie traces"""
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

    """Create Traces for plot"""
    num_trace1 = go.Scatter(
        x=num_mrna["Time"],
        y=num_mrna["mRNA"],
        name="Numeric # of mRNA",
        marker=dict(
            color='rgb(34,163,192)'
        )
    )

    num_trace2 = go.Scatter(
        x=num_protein["Time"],
        y=num_protein['Proteins'],
        name="Numeric # of Protein"
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

    """Plot Numerical vs. ODE sim"""
    ode_num_fig = make_subplots(specs=[[{"secondary_y": True}]])
    # Numerical traces
    ode_num_fig.add_trace(num_trace1, secondary_y=True)
    ode_num_fig.add_trace(num_trace2, secondary_y=False)

    # ODE traces
    ode_num_fig.add_trace(ode_trace1, secondary_y=True,)
    ode_num_fig.add_trace(ode_trace2, secondary_y=False)
    ode_num_fig.update_layout(
        title="ODE and numerical comparison of mRNA and Protein molecules over time",
        xaxis_title="Time (s)",
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        barmode="group",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    ode_num_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", secondary_y=True)
    ode_num_fig.show()
    # ode_num_fig.write_html(num_ode_path, include_plotlyjs=True)
    # ode_num_fig.write_image(num_ode_image_path)

    """Plot Gillespie Data vs ODE Data"""
    gill_fig = make_subplots(specs=[[{"secondary_y": True}]])
    gill_fig.add_trace(gill_trace1, secondary_y=True,)
    gill_fig.add_trace(gill_trace2, secondary_y=False,)
    gill_fig.add_trace(ode_trace1, secondary_y=True,)
    gill_fig.add_trace(ode_trace2, secondary_y=False, )
    gill_fig.update_layout(
        title="Gillespie comparison of mRNA and Protein molecules over time",
        xaxis_title="Time (s)",
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    gill_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", secondary_y=True)
    gill_fig.show()
    # gill_fig.write_html(gill_path)
    # gill_fig.write_image(gill_image_path)

    """Get satistics from dataframes For numerical, gillespie and ODE sims"""
    num_prot_mean, num_prot_var = num_protein.mean(), num_protein.var()
    num_mrna_mean, num_mrna_var = num_mrna.mean(), num_mrna.var()

    ode_prot_mean, ode_prot_var = ode_proteins.mean(), ode_proteins.var()
    ode_mrna_mean, ode_mrna_var = ode_mrna.mean(), ode_mrna.var()

    gill_prot_mean, gill_prot_var = gill_protein.mean(), gill_protein.var()
    gill_mrna_mean, gill_mrna_var = gill_mrna.mean(), gill_mrna.var()

    mrna_mean_df = pandas.DataFrame([
        ["mRNA mean", num_mrna_mean["mRNA"], "Numerical"],
        ["mRNA mean", ode_mrna_mean["mRNA"], "ODE"],
        ["mRNA mean", gill_mrna_mean["mRNA"], "Gillespie"],
    ], columns=["Stat Type", "Value", "Sim Type"])

    mrna_var_df = pandas.DataFrame([
        ["mRNA var", num_mrna_var["mRNA"], "Numerical"],
        ["mRNA var", ode_mrna_var["mRNA"], "ODE"],
        ["mRNA var", gill_mrna_var["mRNA"], "Gillespie"]
    ], columns=["Stat Type", "Value", "Sim Type"])

    prot_mean_df = pandas.DataFrame([
        ["Protein mean", num_prot_mean["Proteins"], "Numerical"],
        ["Protein mean", ode_prot_mean["Proteins"], "ODE"],
        ["Protein mean", gill_prot_mean["Proteins"], "Gillespie"],
    ], columns=["Stat Type", "Value", "Sim Type"])

    prot_var_df = pandas.DataFrame([
        ["Protein var", num_prot_var["Proteins"], "Numerical"],
        ["Protein var", ode_prot_var["Proteins"], "ODE"],
        ["Protein var", gill_prot_var["Proteins"], "Gillespie"],
    ], columns=["Stat Type", "Value", "Sim Type"])

    # stat_fig = make_subplots(rows=2, cols=2)
    mrna_mean_fig = px.bar(mrna_mean_df,
                           x="Stat Type",
                           y="Value",
                           color="Sim Type",
                           barmode="group",
                           labels={
                               "Value": "Average number of Molecules"
                           },
                           title="mRNA mean of different simulations")

    mrna_var_fig = px.bar(mrna_var_df,
                          x="Stat Type",
                          y="Value",
                          color="Sim Type",
                          barmode="group",
                          labels={
                              "Value": "Variance"
                          },
                          title="mRNA variance of different simulations")

    prot_mean_fig = px.bar(prot_mean_df,
                           x="Stat Type",
                           y="Value",
                           color="Sim Type",
                           barmode="group",
                           labels={
                                "Value": "Average number of molecules"
                           },
                           title="Protein Mean of different simulations")

    prot_var_fig = px.bar(prot_var_df,
                          x="Stat Type",
                          y="Value",
                          color="Sim Type",
                          barmode="group",
                          labels={
                              "Value": "Variance"
                          },
                          title="Protein Variance of different simulations",
                          )
    mrna_var_fig.show()
    mrna_mean_fig.show()
    prot_var_fig.show()
    prot_mean_fig.show()

    # stat_ode_num_fig.write_html(stat_path)
    # stat_ode_num_fig.write_image(stat_image_path)


if __name__ == '__main__':
    main()
