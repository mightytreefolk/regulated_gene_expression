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
    numerical_model = UnregulatedGeneExpression(m0=z0[0], p0=z0[1], const=const)

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
    const = [0.0167, 0.167, 0.0022, 0.00125]
    # m0, p0
    initial_conditions = [7.59, 1014.145]

    """Extract data"""
    num_mrna, num_protein = numerical_sim(tmax=tmax, z0=initial_conditions, n=n, const=const)
    ode_mrna, ode_proteins = ode_sim(tmax=tmax, z0=initial_conditions, n=n, const=const)
    gill_model = GillespieUnregulatedGeneExpression(tmax=tmax, m0=initial_conditions[0], p0=initial_conditions[1],
                                                    const=const)
    gill_mrna, gill_protein = gill_model.run_sim()

    """Define plot paths"""
    base_path = "plots/"
    if initial_conditions[0] == 0 and initial_conditions[1] == 0:
        run_name = "from_zero"
    elif initial_conditions[0] == 7.59 and initial_conditions[1] == 6072.675:
        run_name = "steady_state"
    else:
        run_name = "something_different"

    timestamp = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    num_ode_path = os.path.join(base_path, "html/ode_num_compare_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                                run_name=run_name))
    gill_path = os.path.join(base_path, "html/gill_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                  run_name=run_name))
    stat_path = os.path.join(base_path, "html/stat_plot-{run_name}-{time}.html".format(time=timestamp,
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
        name="GILL - mRNA",
        line=dict(color='royalblue', )

    )

    gill_trace2 = go.Scatter(
        x=gill_protein["Time"],
        y=gill_protein['Proteins'],
        name="GILL - Protein",
        line=dict(color='firebrick', )
    )

    """Create Traces for plot"""
    analytical_trace1 = go.Scatter(
        x=num_mrna["Time"],
        y=num_mrna["mRNA"],
        name="Analytical - mRNA",
        line=dict(color='royalblue',)
    )

    analytical_trace2 = go.Scatter(
        x=num_protein["Time"],
        y=num_protein['Proteins'],
        name="Analytical - Protein",
        line=dict(color='firebrick', )

    )

    ode_trace1 = go.Scatter(
        x=ode_mrna["Time"],
        y=ode_mrna["mRNA"],
        name="Numerical - mRNA",
        line=dict(color='royalblue',
                  dash='dash')
    )

    ode_trace2 = go.Scatter(
        x=ode_proteins["Time"],
        y=ode_proteins['Proteins'],
        name="Numerical - Protein",
        line=dict(color='firebrick',
                  dash='dash')

    )

    """Plot Numerical vs. Analytical sim"""
    ode_num_fig = make_subplots(specs=[[{"secondary_y": True}]])
    # Numerical traces
    ode_num_fig.add_trace(analytical_trace1, secondary_y=True)
    ode_num_fig.add_trace(analytical_trace2, secondary_y=False)

    # ODE traces
    ode_num_fig.add_trace(ode_trace1, secondary_y=True,)
    ode_num_fig.add_trace(ode_trace2, secondary_y=False)
    ode_num_fig.update_layout(
        title="Analytical and numerical comparison of mRNA and Protein molecules over time",
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
    ode_num_fig.write_html(num_ode_path, include_plotlyjs=True)
    ode_num_fig.write_image(num_ode_image_path)

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
    gill_fig.write_html(gill_path)
    gill_fig.write_image(gill_image_path)

    # """Get satistics from dataframes For numerical, gillespie and ODE sims"""
    # num_prot_mean, num_prot_var = num_protein.mean(), num_protein.var()
    # num_mrna_mean, num_mrna_var = num_mrna.mean(), num_mrna.var()
    #
    # ode_prot_mean, ode_prot_var = ode_proteins.mean(), ode_proteins.var()
    # ode_mrna_mean, ode_mrna_var = ode_mrna.mean(), ode_mrna.var()
    #
    # gill_prot_mean, gill_prot_var = gill_protein.mean(), gill_protein.var()
    # gill_mrna_mean, gill_mrna_var = gill_mrna.mean(), gill_mrna.var()
    #
    # """Create Stat traces from data"""
    # ode_prot_mean_trace = go.Bar(x=["Protein Mean"],
    #                              y=[ode_prot_mean["Proteins"]],
    #                              name="ODE Protein",
    #                              marker=dict(color=["crimson"])
    #                              )
    # ode_mrna_mean_trace = go.Bar(x=["mRNA Mean"],
    #                              y=[ode_mrna_mean["mRNA"]],
    #                              name="ODE mRNA",
    #                              marker=dict(color=["crimson"])
    #                              )
    # ode_prot_var_trace = go.Bar(x=["Protein Variance"],
    #                             y=[ode_prot_var["Proteins"]],
    #                             name="ODE Protein",
    #                             marker=dict(color=["crimson"])
    #                             )
    # ode_mrna_var_trace = go.Bar(x=["mRNA Variance"],
    #                             y=[ode_mrna_var["mRNA"]],
    #                             name="ODE Protein",
    #                             marker=dict(color=["crimson"])
    #                             )
    #
    # gill_prot_mean_trace = go.Bar(x=["Protein Mean"],
    #                               y=[gill_prot_mean["Proteins"]],
    #                               name="Gillespie Protein",
    #                               marker=dict(color=["orange"])
    #                               )
    # gill_prot_var_trace = go.Bar(x=["Protein Variance"],
    #                              y=[gill_prot_var["Proteins"]],
    #                              name="Gillespie Protein",
    #                              marker=dict(color=["orange"])
    #                              )
    # gill_mrna_mean_trace = go.Bar(x=["mRNA Mean"],
    #                               y=[gill_mrna_mean["mRNA"]],
    #                               name="Gillespie mRNA",
    #                               marker=dict(color=["orange"])
    #                               )
    # gill_mrna_var_trace = go.Bar(x=["mRNA Variance"],
    #                              y=[gill_mrna_var["mRNA"]],
    #                              name="Gillespie mRNA",
    #                              marker=dict(color=["orange"])
    #                              )
    #
    # num_mrna_mean_trace = go.Bar(x=["mRNA Mean"],
    #                              y=[num_mrna_mean["mRNA"]],
    #                              name="Numerical mRNA",
    #                              marker=dict(color=["blue"])
    #                              )
    # num_prot_mean_trace = go.Bar(x=["Protein Mean"],
    #                              y=[num_prot_mean["Proteins"]],
    #                              name="Numerical Protein",
    #                              marker=dict(color=["blue"])
    #                              )
    # num_prot_var_trace = go.Bar(x=["Protein Variance"],
    #                             y=[num_prot_var["Proteins"]],
    #                             name="Numerical Protein",
    #                             marker=dict(color=["blue"])
    #                             )
    # num_mrna_var_trace = go.Bar(x=["mRNA Variance"],
    #                             y=[num_mrna_var["mRNA"]],
    #                             name="Numerical mRNA",
    #                             marker=dict(color=["blue"])
    #                             )
    #
    # """Graph the Stats in a bar chart"""
    # stat_ode_num_fig = make_subplots(rows=2, cols=1,
    #                                  specs=[[{"secondary_y": True}],
    #                                         [{"secondary_y": True}]])
    # stat_ode_num_fig.add_trace(num_prot_mean_trace, row=1, col=1, secondary_y=False)
    # stat_ode_num_fig.add_trace(ode_prot_mean_trace, row=1, col=1, secondary_y=False)
    # stat_ode_num_fig.add_trace(num_mrna_mean_trace, row=1, col=1, secondary_y=True)
    # stat_ode_num_fig.add_trace(ode_mrna_mean_trace, row=1, col=1, secondary_y=True)
    # stat_ode_num_fig.add_trace(gill_prot_mean_trace, row=1, col=1, secondary_y=False)
    # stat_ode_num_fig.add_trace(gill_mrna_mean_trace, row=1, col=1, secondary_y=True)
    #
    # stat_ode_num_fig.add_trace(num_prot_var_trace, row=2, col=1, secondary_y=False)
    # stat_ode_num_fig.add_trace(ode_prot_var_trace, row=2, col=1, secondary_y=False)
    # stat_ode_num_fig.add_trace(num_mrna_var_trace, row=2, col=1, secondary_y=True)
    # stat_ode_num_fig.add_trace(ode_mrna_var_trace, row=2, col=1, secondary_y=True)
    # stat_ode_num_fig.add_trace(gill_prot_var_trace, row=2, col=1, secondary_y=False)
    # stat_ode_num_fig.add_trace(gill_mrna_var_trace, row=2, col=1, secondary_y=True)
    # stat_ode_num_fig.update_layout(
    #     title="Mean and Variance comparisons between numerical, Gillespie and ODE simulations",
    #     yaxis_title="Number of Molecules",
    #     legend_title="Legend",
    #     font=dict(
    #         family="Courier New, monospace",
    #         size=12,
    #         color="Black"
    #     )
    #
    # )
    # stat_ode_num_fig.show()
    # stat_ode_num_fig.write_html(stat_path)
    # stat_ode_num_fig.write_image(stat_image_path)


if __name__ == '__main__':
    main()
