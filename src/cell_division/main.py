import os
import pandas
import plotly.graph_objects as go
import numpy
from plotly.subplots import make_subplots
import glob
from models import Gillespie, CellDivision, DeterministicCellDivision

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


def combine_cell_cycles(sim):
    path = os.path.join(sim, "*.csv")
    for fname in glob.glob(path):
        df = pandas.read_csv(fname, sep='\t')
    result = zip(df["Counter"], df["mRNA"])
    cell_cycles = pandas.DataFrame()
    mrna = []
    l = []
    for i in result:
        if i[0] == 0:
            mrna.append(len(l))
            l.clear()
        else:
            l.append(i[1])
    mrna.pop(0)
    print(mrna)
    print(numpy.array(mrna).mean())




def analytical_plot(run, save):
    """Plot Analytical model"""
    deterministic_run = run.analytical_simulation()
    deterministic_run = deterministic_run.iloc[12600:]
    deterministic_fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02)
    mrna_trace = go.Scatter(x=deterministic_run["Time"],
                            y=deterministic_run["mRNA"],
                            name="mRNA")

    prot_trace = go.Scatter(x=deterministic_run["Time"],
                            y=deterministic_run["Proteins"],
                            name="Proteins")

    deterministic_fig.add_trace(mrna_trace, row=2, col=1)
    deterministic_fig.add_trace(prot_trace, row=1, col=1)
    deterministic_fig.update_layout(
        title="Analytical Cell division comparison of mRNA and Protein molecules over time",
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        barmode="group",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    deterministic_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", row=2, col=1)
    deterministic_fig.update_xaxes(title_text="Time (Hours)", row=2, col=1)
    deterministic_fig.show()
    if save:
        deterministic_fig.write_html("Deterministic.html")
        deterministic_fig.write_image("Deterministic.png")
    else:
        pass


def numerical_plot(run, save):
    """Plot Analytical model"""
    numerical_run = run.numerical_sim()
    numerical_fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02)
    mrna_trace = go.Scatter(x=numerical_run["Time"],
                            y=numerical_run["mRNA"],
                            name="mRNA")

    prot_trace = go.Scatter(x=numerical_run["Time"],
                            y=numerical_run["Proteins"],
                            name="Proteins")

    numerical_fig.add_trace(mrna_trace, row=2, col=1)
    numerical_fig.add_trace(prot_trace, row=1, col=1)
    numerical_fig.update_layout(
        title="Numerical Cell division comparison of mRNA and Protein molecules over time",
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        barmode="group",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    numerical_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", row=2, col=1)
    numerical_fig.update_xaxes(title_text="Time (Hours)", row=2, col=1)
    numerical_fig.show()
    if save:
        numerical_fig.write_html("numerical.html")
        numerical_fig.write_image("numerical.png")
    else:
        pass


def histogram_plot(number_of_runs, sim, save):
    norm_hist = make_subplots(rows=2, cols=1, vertical_spacing=0.02)
    norm_hist.update_layout(
            title="Probability distribution of mRNA and Protein of {n} cell(s)".format(n=number_of_runs),
            yaxis_title="Probability of <b>mRNA</b> Molecules",
            legend_title="Legend",
            font=dict(
                family="Courier New, monospace",
                size=12,
                color="Black"))
    norm_hist.update_yaxes(title_text="Probability of <b>Protein</b>", row=2, col=1)
    norm_hist.update_xaxes(title_text="Number of <b>Molecules</b>", row=2, col=1)
    """Get data from runs"""
    path = os.path.join(sim, "*.csv")
    mrna = []
    prot = []
    for fname in glob.glob(path):
        df = pandas.read_csv(fname, sep='\t')
        mrna.extend(df["mRNA"].tolist())
        prot.extend(df["Proteins"].tolist())
    mrna_hist = go.Histogram(x=mrna, histnorm='probability', name="mRNA Histogram")
    prot_hist = go.Histogram(x=prot, histnorm='probability', name="Protein Histogram")
    norm_hist.add_trace(mrna_hist, row=1, col=1)
    norm_hist.add_trace(prot_hist, row=2, col=1)
    norm_hist.show()
    if save:
        html_file = os.path.join(sim, "hist.html")
        png_file = os.path.join(sim, "hist.png")
        norm_hist.write_html(html_file)
        norm_hist.write_image(png_file)
    else:
        pass


def plot_statistics(number_of_runs, sim, save):
    stat_fig = make_subplots(rows=2, cols=2, vertical_spacing=0.02)
    stat_fig.update_layout(
        title="Mean and Variance for {n} cell(s) dividing".format(n=number_of_runs),
        yaxis_title="Number of Molecules",
        showlegend=False,
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    path = os.path.join(sim, "*.csv")
    mrna_mean = []
    mrna_var = []
    prot_mean = []
    prot_var = []
    for fname in glob.glob(path):
        df = pandas.read_csv(fname, sep='\t')
        mrna_mean.append(df["mRNA"].mean())
        mrna_var.append(df["mRNA"].var())
        prot_mean.append(df["Proteins"].mean())
        prot_var.append(df["Proteins"].var())
    if number_of_runs > 1:
        prot_mean = [numpy.array(prot_mean).mean()]
        prot_var = [numpy.array(prot_var).var()]
        mrna_mean = [numpy.array(mrna_mean).mean()]
        mrna_var = [numpy.array(mrna_var).var()]
    else:
        pass
    gill_protein_mean = go.Bar(x=["Gillespie Protein Mean"],
                               y=prot_mean,
                               text=prot_mean,
                               name="Gillespie Protein",
                               marker=dict(color=["firebrick"]))

    gill_mrna_mean = go.Bar(x=["Gillespie mRNA Mean"],
                            y=mrna_mean,
                            text=mrna_mean,
                            name="Gillespie mRNA",
                            marker=dict(color=["royalblue"]))

    gill_protein_var = go.Bar(x=["Gillespie Protein Variance"],
                              y=prot_var,
                              text=prot_var,
                              name="Gillespie Protein",
                              marker=dict(color=["firebrick"]))
    gill_mrna_var = go.Bar(x=["Gillespie mRNA Variance"],
                           y=mrna_var,
                           text=mrna_var,
                           name="Gillespie mRNA",
                           marker=dict(color=["royalblue"]))
    stat_fig.add_trace(gill_protein_mean, row=1, col=1)
    stat_fig.add_trace(gill_mrna_mean, row=1, col=2)
    stat_fig.add_trace(gill_protein_var, row=2, col=1)
    stat_fig.add_trace(gill_mrna_var, row=2, col=2)
    stat_fig.show()
    if save:
        html_file = os.path.join(sim, "stats.html")
        png_file = os.path.join(sim, "stats.png")
        stat_fig.write_html(html_file)
        stat_fig.write_image(png_file)
    else:
        pass


def plot_gillespie(number_of_runs, sim, save):
    cell_div_fig = make_subplots(specs=[[{"secondary_y": True}],
                                        [{"secondary_y": False}]],
                                 rows=2,
                                 cols=1,
                                 row_heights=[0.8, 0.2],
                                 shared_xaxes=True,
                                 vertical_spacing=0.02)
    cell_div_fig.update_layout(
        title="Cell division comparison of mRNA and Protein molecules over time",
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    cell_div_fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", secondary_y=True)
    cell_div_fig.update_xaxes(title_text="Time (Hours)", row=2, col=1)
    cell_div_fig.update_yaxes(title_text="Number of <b>Genes</b>", row=2, col=1)

    path = os.path.join(sim, "*.csv")
    if number_of_runs == 1:
        for fname in glob.glob(path):
            df = pandas.read_csv(fname, sep='\t')
        """Gillespie model of cell division"""
        mrna_trace = go.Scatter(x=df["Time"],
                                y=df["mRNA"],
                                name="mRNA",
                                line=dict(color='royalblue', )
                                )
        protein_trace = go.Scatter(x=df["Time"],
                                   y=df["Proteins"],
                                   name="Protein",
                                   line=dict(color='firebrick', )
                                   )

        genes_trace = go.Scatter(x=df["Time"],
                                 y=df["Gene Number"],
                                 name="Number of genes")

        cell_div_fig.add_trace(mrna_trace, secondary_y=True, row=1, col=1)
        cell_div_fig.add_trace(protein_trace, secondary_y=False, row=1, col=1)
        cell_div_fig.add_trace(genes_trace, row=2, col=1)
        for i in division(df):
            cell_div_fig.add_shape(i, row=1, col=1, )
        # cell_div_fig.update_shapes(dict(xref='x', yref='y'))
        cell_div_fig.show()
    else:
        pass
    if save:
        html_path = os.path.join(sim, "Gillespie.html")
        image_path = os.path.join(sim, "Gillespie.png")
        cell_div_fig.write_image(image_path)
        cell_div_fig.write_html(html_path)
    else:
        pass

def main():
    """Constants to be changed by user"""
    # seconds for sim
    tmax = 43200
    number_of_datapoints = 43200
    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 0.167, 0.0022, 0]

    # m0, p0 [0, 0]
    initial_conditions = [7, 1014]

    number_of_simulations = 1

    save = True

    # """Initiate Numerical sim"""
    # numerical = DeterministicCellDivision(tmax=tmax,
    #                                       num_of_datapoints=number_of_datapoints,
    #                                       m0=initial_conditions[0],
    #                                       p0=initial_conditions[1],
    #                                       const=const)
    # numerical_plot(numerical, save)
    #
    # """Initiate Analytical sim"""
    # deterministic = DeterministicCellDivision(tmax=tmax,
    #                                           num_of_datapoints=number_of_datapoints,
    #                                           m0=initial_conditions[0],
    #                                           p0=initial_conditions[1],
    #                                           const=const)
    # analytical_plot(deterministic, save=save)

    """Begin Gillespie Simulation"""
    gillespie_cell_model = CellDivision(tmax=tmax, m0=initial_conditions[0], p0=initial_conditions[1], const=const,
                                        number_of_sims=number_of_simulations)

    run = gillespie_cell_model.multiple_cells()
    combine_cell_cycles(sim=run)
    #
    # """Different plots for the Gillespie data"""
    # histogram_plot(number_of_simulations, sim=run, save=save)
    # plot_statistics(number_of_simulations, sim=run, save=save)
    # plot_gillespie(number_of_runs=number_of_simulations, sim=run, save=save)


if __name__ == '__main__':
    main()
