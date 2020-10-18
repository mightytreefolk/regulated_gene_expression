import numpy as np
import os
import math
import pandas
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime
from models import Gillespie, CellDivision

n_A = 6.023E23  # Avogadro's Number
e_coli_vol = 6.5E-16  # Liters


def division(df):
    result = zip(df["Time"], df["Counter"])
    traces = []
    for i in result:
        if i[1] == 0:
            trace = dict(
                type="line",
                x0=i[0],
                y0=0,
                x1=i[0],
                y1=2500,
                line=dict(
                    color="Black",
                    width=1,
                    dash="dashdot"
                ))
            traces.append(trace)
        else:
            pass
    return traces



def main():
    # seconds for sim
    tmax = 43200

    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 0.167, 0.0022, 0.00125]

    # m0, p0 [0, 0]
    initial_conditions = [7, 1014]

    gillespie_cell_model = CellDivision(tmax=tmax, m0=initial_conditions[0], p0=initial_conditions[1], const=const)
    run = gillespie_cell_model.sim()
    run = run.iloc[40000:]

    mrna_trace = go.Scatter(x=run["Time"],
                            y=run["mRNA"],
                            name="mRNA",
                            line=dict(color='royalblue',)
                            )
    protein_trace = go.Scatter(x=run["Time"],
                               y=run["Proteins"],
                               name="Protein",
                               line=dict(color='firebrick', )
                               )

    genes_trace = go.Scatter(x=run["Time"],
                             y=run["Gene Number"],
                             name="Number of genes")

    cell_div_fig = make_subplots(specs=[[{"secondary_y": True}]])
    cell_div_fig.add_trace(mrna_trace, secondary_y=True)
    cell_div_fig.add_trace(protein_trace, secondary_y=False)
    cell_div_fig.add_trace(genes_trace, secondary_y=True)
    for i in division(run):
        cell_div_fig.add_shape(i, secondary_y=True)

    cell_div_fig.update_layout(
        title="Cell division comparison of mRNA and Protein molecules over time",
        xaxis_title="Time (Hours)",
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        barmode="group",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    cell_div_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", secondary_y=True)
    cell_div_fig.update_shapes(dict(xref='x', yref='y'))
    cell_div_fig.show()

    """Plot Histogram of mRNA"""
    norm_mrna_hist = go.Figure()
    hist = go.Histogram(x=run["mRNA"], histnorm='probability', name="mRNA Histogram")
    norm_mrna_hist.add_trace(hist)
    norm_mrna_hist.update_layout(
        title="Probability distribution of mRNA",
        xaxis_title="Number of Molecules",
        yaxis_title="Probability of <b>mRNA</b> Molecules",
        legend_title="Legend",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    norm_mrna_hist.show()

    """Plot statistics"""
    gill_protein_mean = go.Bar(x=["Gillespie Protein Mean"],
                               y=[run["Proteins"].mean()],
                               name="Gillespie Protein",
                               marker=dict(color=["firebrick"]))

    gill_mrna_mean = go.Bar(x=["Gillespie mRNA Mean"],
                            y=[run["mRNA"].mean()],
                            name="Gillespie mRNA",
                            marker=dict(color=["royalblue"]))

    gill_protein_var = go.Bar(x=["Gillespie Protein Variance"],
                              y=[run["Proteins"].var()],
                              name="Gillespie Protein",
                              marker=dict(color=["firebrick"]))
    gill_mrna_var = go.Bar(x=["Gillespie mRNA Variance"],
                           y=[run["mRNA"].var()],
                           name="Gillespie mRNA",
                           marker=dict(color=["royalblue"]))

    stat_fig = make_subplots(rows=2, cols=2)
    stat_fig.add_trace(gill_protein_mean, row=1, col=1)
    stat_fig.add_trace(gill_mrna_mean, row=1, col=2)
    stat_fig.add_trace(gill_protein_var, row=2, col=1)
    stat_fig.add_trace(gill_mrna_var, row=2, col=2)
    stat_fig.update_layout(
        title="Mean and Variance for cells dividing",
        yaxis_title="Number of Molecules",
        showlegend=False,
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    stat_fig.show()

if __name__ == '__main__':
    main()