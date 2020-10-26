import os
import pandas
import uuid
import math
import plotly.graph_objects as go
from scipy.stats import sem
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
                y1=25,
                line=dict(
                    color="Black",
                    width=1,
                    dash="dashdot"
                ))
            traces.append(trace)
        else:
            pass
    return traces


def gaussian(u, s, x):
    gauss = (1/(s * math.sqrt(2 * math.pi))) * math.exp(-0.5 * ((x-u)/s)**2)
    return gauss


def average_cycle_times(sim, save):
    path = os.path.join(sim, "*.csv")
    for fname in glob.glob(path):
        df = pandas.read_csv(fname, sep='\t')
    time_at_division = []
    time_at_two_genes = []
    result = zip(df["Clock"].tolist(), df["Counter"].tolist())
    for i in result:
        if i[1] == 3600:
            time_at_division.append(i[0])
        elif i[1] == 1800:
            time_at_two_genes.append(i[0])

    time_at_division = [x / 60 for x in time_at_division]
    time_at_two_genes =[x / 60 for x in time_at_two_genes]
    print("Max time division: ", max(time_at_division))
    print("Min time division: ", min(time_at_division))
    print("Division mean: ", numpy.array(time_at_division).mean())
    print("Division std: ", numpy.array(time_at_division).std())

    print("Max time at 2 genes: ", max(time_at_two_genes))
    print("Min time at 2 genes: ", min(time_at_two_genes))
    print("Two genes mean: ", numpy.array(time_at_two_genes).mean())
    print("Two genes std: ", numpy.array(time_at_two_genes).std())
    genes_mean = numpy.array(time_at_two_genes).mean()
    genes_std = numpy.array(time_at_two_genes).std()
    divide_mean = numpy.array(time_at_division).mean()
    divide_std = numpy.array(time_at_division).std()
    gauss_divide = []
    gauss_genes = []
    for i in numpy.sort(time_at_division):
        gauss_divide.append(gaussian(u=divide_mean, s=divide_std, x=i))
    for i in numpy.sort(time_at_two_genes):
        gauss_genes.append(gaussian(u=genes_mean, s=genes_std, x=i))

    fig = make_subplots(rows=2, cols=1)
    fig.update_yaxes(title_text="Probability of <b>Two Genes</b> (minutes)", row=2, col=1)
    fig.update_xaxes(title_text="Time (Min)")
    fig.update_layout(
        title="Probability of time of division and two genes in cell ({n} cell cycles)".format(n=len(time_at_division)),
        yaxis_title="Probability of <b>division</b> (minutes)",
        legend_title="Legend",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"))

    division_trace = go.Histogram(x=time_at_division, histnorm='probability', name="Division Histogram")
    division_gauss_trace = go.Scatter(x=numpy.sort(time_at_division), y=gauss_divide,  showlegend=False)
    two_gene_trace = go.Histogram(x=time_at_two_genes, histnorm='probability', name="2 Genes Histogram")
    two_gene_gauss_trace = go.Scatter(x=numpy.sort(time_at_two_genes), y=gauss_genes,  showlegend=False)
    fig.add_trace(division_trace, row=1, col=1)
    fig.add_trace(division_gauss_trace, row=1, col=1,)
    fig.add_trace(two_gene_trace, row=2, col=1)
    fig.add_trace(two_gene_gauss_trace, row=2, col=1)
    fig.show()
    if save:
        html = os.path.join(sim, "Average_times.html")
        image = os.path.join(sim, "Average_times.png")
        fig.write_html(html)
        fig.write_image(image)
    else:
        pass


def combine_cell_cycles(sim, save, const):
    path = os.path.join(sim, "*.csv")
    for fname in glob.glob(path):
        df = pandas.read_csv(fname, sep='\t')
    # df = df.iloc[25200:]
    ml = []
    tl = []
    pl = []
    number_of_cycles = []
    mrna_cell_cycles = pandas.DataFrame()
    protein_cell_cycles = pandas.DataFrame()
    counter = pandas.DataFrame()
    result = zip(df["Counter"].tolist(), df["mRNA"].tolist(), df["Divide"].tolist(), df["Proteins"].tolist())
    for i in result:
        if i[2] == "Yes":
            x = str(uuid.uuid4())
            mrna_cell_cycles["Cycle_{}".format(x)] = ml
            protein_cell_cycles["Cycle_{}".format(x)] = pl
            counter["Cycle_{}".format(x)] = tl
            number_of_cycles.append(i[0])
            ml.clear()
            tl.clear()
            pl.clear()
            ml.append(i[1])
            tl.append(i[0])
            pl.append(i[3])
        else:
            ml.append(i[1])
            tl.append(i[0])
            pl.append(i[3])

    mrna_cell_cycles["Average"] = mrna_cell_cycles.mean(axis=1)
    mrna_cell_cycles["std"] = mrna_cell_cycles.std(axis=1)
    mrna_cell_cycles["Counter"] = numpy.linspace(0, 1, len(mrna_cell_cycles["Average"].tolist()))
    protein_cell_cycles["Average"] = protein_cell_cycles.mean(axis=1)
    protein_cell_cycles["std"] = protein_cell_cycles.std(axis=1)
    protein_cell_cycles["Counter"] = numpy.linspace(0, 1, len(protein_cell_cycles["Average"].tolist()))

    mrna_std_minus = []
    mrna_std_plus = []

    prot_std_minus = []
    prot_std_plus = []

    mrna_stats = zip(mrna_cell_cycles["Average"].tolist(), mrna_cell_cycles["std"].tolist())
    prot_stats = zip(protein_cell_cycles["Average"].tolist(), protein_cell_cycles["std"].tolist())

    for i, j in mrna_stats:
        mrna_std_minus.append(i-j)
        mrna_std_plus.append(i+j)
    for i, j in prot_stats:
        prot_std_minus.append(i-j)
        prot_std_plus.append(i+j)

    numerical = DeterministicCellDivision(tmax=3600,
                                          num_of_datapoints=3600,
                                          m0=7.59,
                                          p0=1014.145,
                                          const=const)
    numerical_run = numerical.numerical_sim()

    """Just graphing things"""
    fig = make_subplots(rows=2, cols=1, vertical_spacing=0.02, shared_xaxes=True)
    # Making traces for plot
    mrna_std_minus_trace = go.Scatter(x=mrna_cell_cycles["Counter"],
                                      y=mrna_std_minus,
                                      name="Gillespie mRNA STD",
                                      line=dict(color='darkgrey'),
                                      fill=None,
                                      showlegend=False)
    mrna_std_plus_trace = go.Scatter(x=mrna_cell_cycles["Counter"],
                                     y=mrna_std_plus,
                                     name="Gillespie mRNA STD",
                                     line=dict(color='darkgrey'),
                                     fill='tonexty',
                                     showlegend=False)

    mrna_trace = go.Scatter(x=mrna_cell_cycles["Counter"],
                            y=mrna_cell_cycles["Average"],
                            name="Gillespie mRNA",
                            line=dict(color='royalblue'))
    numerical_mrna_trace = go.Scatter(x=numerical_run["Time"],
                                      y=numerical_run["mRNA"],
                                      name="Numerical mRNA",
                                      line=dict(color='royalblue',
                                                dash='dash'))
    prot_std_minus_trace = go.Scatter(x=protein_cell_cycles["Counter"],
                                      y=prot_std_minus,
                                      name="Gillespie Protein STD",
                                      line=dict(color='darkgrey'),
                                      fill=None,
                                      showlegend=False)
    prot_std_plus_trace = go.Scatter(x=protein_cell_cycles["Counter"],
                                     y=prot_std_plus,
                                     name="Gillespie Protein STD",
                                     line=dict(color='darkgrey'),
                                     fill='tonexty',
                                     showlegend=False)

    prot_trace = go.Scatter(x=protein_cell_cycles["Counter"],
                            y=protein_cell_cycles["Average"],
                            name="Gillespie Protein",
                            line=dict(color='firebrick'))
    numerical_prot_trace = go.Scatter(x=numerical_run["Time"],
                                      y=numerical_run["Proteins"],
                                      name="Numerical Proteins",
                                      line=dict(color='firebrick',
                                                dash='dash'))
    fig.add_trace(prot_trace, row=1, col=1)
    fig.add_trace(numerical_prot_trace, row=1, col=1)
    fig.add_trace(prot_std_minus_trace, row=1, col=1)
    fig.add_trace(prot_std_plus_trace, row=1, col=1,)

    fig.add_trace(numerical_mrna_trace, row=2, col=1)
    fig.add_trace(mrna_trace, row=2, col=1)
    fig.add_trace(mrna_std_minus_trace, row=2, col=1)
    fig.add_trace(mrna_std_plus_trace, row=2, col=1)
    fig.update_yaxes(title_text="Number of <b>mRNA</b> Molecules", row=2, col=1)
    fig.update_xaxes(title_text="Time (Hours)", row=2, col=1)
    fig.update_layout(
        title="Cell cycle average of Protein and mRNA over time for {n} cell cycles".format(n=len(number_of_cycles)),
        yaxis_title="Number of <b>Protein</b> Molecules",
        legend_title="Legend",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    fig.show()
    if save:
        html = os.path.join(sim, "Cell_cycle_average.html")
        image = os.path.join(sim, "Cell_cycle_average.png")
        fig.write_html(html)
        fig.write_image(image)
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


def histogram_plot(sim, save):
    norm_hist = make_subplots(rows=2, cols=1, vertical_spacing=0.02)
    norm_hist.update_layout(
            title="Probability distribution of mRNA and Protein",
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


def plot_statistics(sim, save, const):
    stat_fig = make_subplots(rows=2, cols=2, vertical_spacing=0.02)
    path = os.path.join(sim, "*.csv")
    mrna_cell_cycles = pandas.DataFrame()
    protein_cell_cycles = pandas.DataFrame()
    ml = []
    pl = []
    for fname in glob.glob(path):
        df = pandas.read_csv(fname, sep='\t')
        result = zip(df["Counter"].tolist(), df["mRNA"].tolist(), df["Divide"].tolist(), df["Proteins"].tolist())
        for i in result:
            if i[2] == "Yes":
                x = str(uuid.uuid4())
                mrna_cell_cycles["Cycle_{}".format(x)] = ml
                protein_cell_cycles["Cycle_{}".format(x)] = pl
                ml.clear()
                pl.clear()
                ml.append(i[1])
                pl.append(i[3])
            else:
                ml.append(i[1])
                pl.append(i[3])
    mrna_mean = mrna_cell_cycles.mean(axis=0).tolist()
    mrna_var = mrna_cell_cycles.var(axis=0).tolist()
    prot_mean = protein_cell_cycles.mean(axis=0).tolist()
    prot_var = protein_cell_cycles.var(axis=0).tolist()
    mrna_sem = sem(mrna_mean)
    prot_sem = sem(prot_mean)

    numerical = DeterministicCellDivision(tmax=3600,
                                          num_of_datapoints=3600,
                                          m0=7.59,
                                          p0=1014.145,
                                          const=const)
    numerical_run = numerical.numerical_sim()
    num_protein_mean = go.Bar(x=["Numerical Protein Mean"],
                              y=[numerical_run["Proteins"].mean()],
                              name="Numerical Protein",
                              marker=dict(color=["darkgrey"]))

    num_protein_var = go.Bar(x=["Numerical Protein Variance"],
                             y=[numerical_run["Proteins"].var()],
                             name="Numerical Protein",
                             marker=dict(color=["darkgrey"]))

    num_mrna_mean = go.Bar(x=["Numerical mRNA Mean"],
                           y=[numerical_run["mRNA"].mean()],
                           name="Numerical mRNA",
                           marker=dict(color=["darkgrey"]))
    num_mrna_var = go.Bar(x=["Numerical mRNA Variance"],
                          y=[numerical_run["mRNA"].var()],
                          name="Numerical mRNA",
                          marker=dict(color=["darkgrey"]))

    gill_protein_mean = go.Bar(x=["Gillespie Protein Mean"],
                               y=[numpy.array(prot_mean).mean()],
                               text=prot_mean,
                               name="Gillespie Protein",
                               marker=dict(color=["firebrick"]),
                               error_y=dict(type='data', array=[prot_sem])
                               )

    gill_mrna_mean = go.Bar(x=["Gillespie mRNA Mean"],
                            y=[numpy.array(mrna_mean).mean()],
                            text=mrna_mean,
                            name="Gillespie mRNA",
                            marker=dict(color=["royalblue"]),
                            error_y=dict(type='data', array=[mrna_sem])
                            )

    gill_protein_var = go.Bar(x=["Gillespie Protein Variance"],
                              y=[numpy.array(prot_var).mean()],
                              text=prot_var,
                              name="Gillespie Protein",
                              marker=dict(color=["firebrick"]),
                              error_y=dict(type='data', array=[numpy.array(prot_var).std()]))
    gill_mrna_var = go.Bar(x=["Gillespie mRNA Variance"],
                           y=[numpy.array(mrna_var).mean()],
                           text=mrna_var,
                           name="Gillespie mRNA",
                           marker=dict(color=["royalblue"]),
                           error_y=dict(type='data', array=[numpy.array(mrna_var).std()]))

    stat_fig.add_trace(gill_protein_mean, row=1, col=1)
    stat_fig.add_trace(num_protein_mean, row=1, col=1)

    stat_fig.add_trace(gill_mrna_mean, row=1, col=2)
    stat_fig.add_trace(num_mrna_mean, row=1, col=2)

    stat_fig.add_trace(gill_protein_var, row=2, col=1)
    stat_fig.add_trace(num_protein_var, row=2, col=1)

    stat_fig.add_trace(gill_mrna_var, row=2, col=2)
    stat_fig.add_trace(num_mrna_var, row=2, col=2)
    stat_fig.update_layout(
        title="Mean and Variance for dividing for {n} cell cycles".format(n=len(mrna_mean)),
        yaxis_title="Number of Molecules",
        showlegend=False,
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
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
            cell_div_fig.add_shape(i, row=1, col=1, secondary_y=True)
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
    tmax = 360000
    number_of_datapoints = 43200
    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 0.167, 0.0022, 0.0]

    # m0, p0 [0, 0]
    initial_conditions = [7, 1014]

    number_of_simulations = 1
    # time points (in unites of dt)
    cell_cycle = 6300
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
                                        number_of_sims=number_of_simulations, cell_cycle=cell_cycle)
    run = gillespie_cell_model.multiple_cells()
    combine_cell_cycles(sim=run, save=save, const=const)
    # average_cycle_times(sim=run, save=save)

    # """Different plots for the Gillespie data"""
    # histogram_plot(sim=run, save=save)
    plot_statistics(sim=run, save=save, const=const)
    # plot_gillespie(number_of_runs=number_of_simulations, sim=run, save=save)


if __name__ == '__main__':
    main()
