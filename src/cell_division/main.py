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
    cell_div_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", secondary_y=True)
    cell_div_fig.update_shapes(dict(xref='x', yref='y'))
    cell_div_fig.show()


if __name__ == '__main__':
    main()