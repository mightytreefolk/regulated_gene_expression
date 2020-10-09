import numpy as np
import os
import pandas
import plotly.express as px
import plotly.graph_objects as go
import json
from plotly.subplots import make_subplots
from datetime import datetime
from scipy.integrate import odeint
from model import UnregulatedGeneExpression, GillespieUnregulatedGeneExpression

n_A = 6.023E23  # Avogadro's Number
e_coli_vol = 6.5E-16  # Liters

'''
Arguments are as follows:
Dataframe1 (mRNA), Dataframe2 (Proteins, number of cells
'''


def gillespie_traces(mrna, prot, n, plot_average):
    if n == 1:
        gill_trace1 = go.Scatter(
            x=mrna["Time"],
            y=mrna["mRNA"],
            name="GILL - mRNA",
            line=dict(color='royalblue', ))
        gill_trace2 = go.Scatter(
            x=prot["Time"],
            y=prot['Proteins'],
            name="GILL - Protein",
            line=dict(color='firebrick', ))
        return gill_trace1, gill_trace2
    elif n > 1 and plot_average:
        gill_trace1 = go.Scatter(
            x=mrna["Average_Time"],
            y=mrna["Average"],
            name="GILL - mRNA",
            line=dict(color='royalblue', ))
        gill_trace2 = go.Scatter(
            x=prot["Average_Time"],
            y=prot['Average'],
            name="GILL - Protein",
            line=dict(color='firebrick', ))
        return [gill_trace1, gill_trace2]
    elif n > 1 and not plot_average:
        traces = []
        for i in range(0, n):
            mrna_trace = go.Scatter(
                x=mrna["mRNA_Run_time{t}".format(t=i)],
                y=mrna["Run{n}".format(n=i)],
                name="GILL - mRNA",
                line=dict(color='royalblue', ))
            prot_trace = go.Scatter(
                x=prot["prot_Run_time{t}".format(t=i)],
                y=prot["Run{n}".format(n=i)],
                name="GILL - Protein",
                line=dict(color='firebrick', ))
            traces.append(mrna_trace)
            traces.append(prot_trace)
        return traces


def main():
    # seconds for sims (for all sims)
    tmax = 10000
    # number of data points (For numerical and analytical)
    n = 10000
    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 0.167, 0.0022, 0.00125]
    # m0, p0
    initial_conditions = [7.59, 1014.145]

    number_of_cells = 30

    plot_average = True

    """Extract data"""
    analytical_numerical_model = UnregulatedGeneExpression(tmax=tmax,
                                                           num_of_datapoints=n,
                                                           m0=initial_conditions[0],
                                                           p0=initial_conditions[1],
                                                           const=const)

    analytical_mrna, analytical_prot = analytical_numerical_model.analytical_sim()

    numerical_mrna, numerical_proteins = analytical_numerical_model.numerical_sim()

    gill_model = GillespieUnregulatedGeneExpression(tmax=tmax,
                                                    m0=initial_conditions[0],
                                                    p0=initial_conditions[1],
                                                    const=const,
                                                    num_cells=number_of_cells)

    gill_mrna, gill_protein = gill_model.multiple_cells_sim()
    traces = gillespie_traces(gill_mrna, gill_protein, number_of_cells, plot_average)



    """Define plot paths"""
    base_path = "plots/"
    if initial_conditions[0] == 0 and initial_conditions[1] == 0:
        run_name = "from_zero"
    elif initial_conditions[0] == 7.59 and initial_conditions[1] == 1014.145:
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

    """Create Traces for plot"""
    analytical_trace1 = go.Scatter(
        x=analytical_mrna["Time"],
        y=analytical_mrna["mRNA"],
        name="Analytical - mRNA",
        line=dict(color='royalblue',)
    )

    analytical_trace2 = go.Scatter(
        x=analytical_prot["Time"],
        y=analytical_prot['Proteins'],
        name="Analytical - Protein",
        line=dict(color='firebrick', )

    )

    numerical_trace1 = go.Scatter(
        x=numerical_mrna["Time"],
        y=numerical_mrna["mRNA"],
        name="Numerical - mRNA",
        line=dict(color='royalblue',
                  dash='dash')
    )

    numerical_trace2 = go.Scatter(
        x=numerical_proteins["Time"],
        y=numerical_proteins['Proteins'],
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
    ode_num_fig.add_trace(numerical_trace1, secondary_y=True,)
    ode_num_fig.add_trace(numerical_trace2, secondary_y=False)
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
    # ode_num_fig.write_html(num_ode_path, include_plotlyjs=True)
    # ode_num_fig.write_image(num_ode_image_path)

    """Plot Gillespie Data vs ODE Data"""
    gill_fig = make_subplots(specs=[[{"secondary_y": True}]])
    for i in traces:
        j = i.to_plotly_json()
        if j['name'] == 'GILL - mRNA':
            gill_fig.add_trace(i, secondary_y=True)
        else:
            gill_fig.add_trace(i, secondary_y=False)
    gill_fig.add_trace(numerical_trace1, secondary_y=True,)
    gill_fig.add_trace(numerical_trace2, secondary_y=False, )
    gill_fig.update_layout(
        title="Gillespie comparison of mRNA and Protein molecules over time for {n} cells".format(n=number_of_cells),
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
