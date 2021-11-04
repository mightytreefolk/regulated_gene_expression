import numpy as np
import os
import math
import plotly.graph_objects as go

from scipy.stats import sem
from plotly.subplots import make_subplots
from datetime import datetime
from models import UnregulatedGeneExpression, GillespieUnregulatedGeneExpression

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

"""
Input is as follows:
v is a vector of the number of molecules (i.e. mRNA or Prot)
k is the creation rate
d is the degradation rate
"""
def prob_dist(v, k, d):
    prob_vector = []
    v.sort()
    largest = v[-1]
    for n in range(0, round(largest)):
        p = math.e**(-k/d) * ((k/d)**round(n))/math.factorial(round(n))
        prob_vector.append(p)
    return [list(range(0, round(largest))), prob_vector]


def main():
    # seconds for sims (for all sims)
    tmax = 10000
    # number of data points (For numerical and analytical)
    n = 10000
    # k0 (mRNA), k1 (protein), dm, dp
    const = [0.0167, 0.167, 0.0022, 0.00125]
    # m0, p0 [0, 0]
    initial_conditions = [7, 1014]

    number_of_cells = 10

    plot_average = True
    save = True

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
    num_ode_image_path = os.path.join(base_path,"images/ode_num_compare_plot-{run_name}-{time}.png".format(time=timestamp,
                                                                                                           run_name=run_name))

    gill_path = os.path.join(base_path, "html/gill_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                       run_name=run_name))

    gill_image_path = os.path.join(base_path, "images/gill_plot-{run_name}-{time}.png".format(time=timestamp,
                                                                                              run_name=run_name))

    stat_image_path = os.path.join(base_path, "images/stat_plot-{run_name}-{time}.png".format(time=timestamp,
                                                                                              run_name=run_name))
    stat_path = os.path.join(base_path, "html/stat_plot-{run_name}-{time}.html".format(time=timestamp,
                                                                                       run_name=run_name))

    norm_mrna_hist_path = os.path.join(base_path, "html/Histogram-{run_name}-{time}.html".format(time=timestamp,
                                                                                                 run_name=run_name))

    norm_mrna_hist_image_path = os.path.join(base_path, "images/Histogram-{run_name}-{time}.png".format(time=timestamp,
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
    # Create figure for traces
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

    """Plot Gillespie Data vs ODE Data"""
    gill_fig = make_subplots(specs=[[{"secondary_y": True}]])
    gill_traces = gillespie_traces(gill_mrna, gill_protein, number_of_cells, plot_average)
    for i in gill_traces:
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

    """Create Histogram from mRNA data"""
    if number_of_cells == 1:
        norm_mrna_hist = go.Figure()
        mrna_prob_data = prob_dist(np.array(gill_mrna["mRNA"]), const[0], const[2])
        mrna_dist = go.Scatter(x=mrna_prob_data[0], y=mrna_prob_data[1], name="Probability distribution")
        hist = go.Histogram(x=gill_mrna["mRNA"], histnorm='probability', name="mRNA Histogram")
        norm_mrna_hist.add_trace(mrna_dist)
        norm_mrna_hist.add_trace(hist)
        norm_mrna_hist.update_layout(
            title="Probability distribution of mRNA for {n} cells".format(n=number_of_cells),
            xaxis_title="Number of Molecules",
            yaxis_title="Probability of <b>mRNA</b> Molecules",
            legend_title="Legend",
            font=dict(
                family="Courier New, monospace",
                size=12,
                color="Black"
            )
        )

    else:
        total_mrna = []
        for i in range(0, number_of_cells):
            mrna = gill_mrna["Run{n}".format(n=i)]
            for m in mrna:
                total_mrna.append(m)
        norm_mrna_hist = go.Figure()
        mrna_prob_data = prob_dist(total_mrna, const[0], const[2])
        mrna_dist = go.Scatter(x=mrna_prob_data[0], y=mrna_prob_data[1], name="Probability distribution")
        hist = go.Histogram(x=total_mrna, histnorm='probability density', name="mRNA Histogram")
        norm_mrna_hist.add_trace(mrna_dist)
        norm_mrna_hist.add_trace(hist)
        norm_mrna_hist.update_layout(
            title="Probability distribution of mRNA for {n} cells".format(n=number_of_cells),
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

    """Get satistics from dataframes For gillespie and ODE sims"""
    num_prot_mean, ode_prot_var = numerical_proteins.mean(), numerical_proteins.var()
    num_mrna_mean, ode_mrna_var = numerical_mrna.mean(), numerical_mrna.var()

    if number_of_cells == 1:
        gill_prot_mean, gill_prot_var = gill_protein["Proteins"].mean(), gill_protein["Proteins"].var()
        gill_mrna_mean, gill_mrna_var = gill_mrna["mRNA"].mean(), gill_mrna["mRNA"].var()
    else:
        gill_prot_var = []
        gill_mrna_var = []
        """
        To Calculate the error for the variance, I made a list of the variances from each run. I then took the mean 
        of the list of variances and used that as my variance for the simulation. To get the error, I took the std of 
        the list of variances.
        """
        for i in range(0, number_of_cells):
            gill_prot_var.append(gill_protein["Run{n}".format(n=i)].var())
            gill_mrna_var.append(gill_mrna["Run{n}".format(n=i)].var())
        gill_prot_var = np.array(gill_prot_var)
        gill_mrna_var = np.array(gill_mrna_var)
        gill_prot_mean= gill_protein["Average"].mean()
        gill_mrna_mean= gill_mrna["Average"].mean()
        prot_sem = sem(gill_protein["Average"].tolist())
        mrna_sem = sem(gill_mrna["Average"].tolist())

    """Create Stat traces from data"""
    num_prot_mean_trace = go.Bar(x=["Numerical Protein Mean"],
                                 y=[num_prot_mean["Proteins"]],
                                 name="ODE Protein",
                                 marker=dict(color=["crimson"])
                                 )
    num_mrna_mean_trace = go.Bar(x=["Numerical mRNA Mean"],
                                 y=[num_mrna_mean["mRNA"]],
                                 name="ODE mRNA",
                                 marker=dict(color=["crimson"])
                                 )

    gill_prot_mean_trace = go.Bar(x=["Gillespie Protein Mean"],
                                  y=[gill_prot_mean],
                                  name="Gillespie Protein",
                                  marker=dict(color=["orange"]),
                                  error_y=dict(type='data', array=[prot_sem])
                                  )
    gill_prot_var_trace = go.Bar(x=["Gillespie Protein Variance"],
                                 y=[gill_prot_var.mean()],
                                 name="Gillespie Protein",
                                 marker=dict(color=["orange"]),
                                 error_y=dict(type='data', array=[sem(gill_prot_var)])
                                 )
    gill_mrna_mean_trace = go.Bar(x=["Gillespie mRNA Mean"],
                                  y=[gill_mrna_mean],
                                  name="Gillespie mRNA",
                                  marker=dict(color=["orange"]),
                                  error_y=dict(type='data', array=[mrna_sem])
                                  )
    gill_mrna_var_trace = go.Bar(x=["Gillespie mRNA Variance"],
                                 y=[gill_mrna_var.mean()],
                                 name="Gillespie mRNA",
                                 marker=dict(color=["orange"]),
                                 error_y=dict(type='data', array=[sem(gill_mrna_var)])
                                 )
    theoretical_mrna_mean = go.Bar(x=["Theoretical mRNA Mean"],
                                   y=[const[0]/const[2]],
                                   name="Theoretical mRNA Mean",
                                   marker=dict(color=["darkgrey"])
                                   )

    theoretical_mrna_var = go.Bar(x=["Theoretical mRNA Variance"],
                                  y=[const[0]/const[2]],
                                  name="Theoretical mRNA Variance",
                                  marker=dict(color=["darkgrey"])
                                  )
    # theoretical_prot_var = go.Bar(x=["Theoretical Protein Variance"],
    #                               y=[(const[0]*const[1])/(const[2]*const[3])],
    #                               name="Theoretical Protein Variance",
    #                               marker=dict(color=["darkgrey"])
    #                               )
    theoretical_prot_mean = go.Bar(x=["Theoretical Protein Mean"],
                                  y=[(const[0]*const[1])/(const[2]*const[3])],
                                  name="Theoretical Protein Mean",
                                  marker=dict(color=["darkgrey"])
                                  )


    """Graph the Stats in a bar chart"""
    stat_ode_num_fig = make_subplots(rows=2, cols=2)
    stat_ode_num_fig.add_trace(num_prot_mean_trace, row=1, col=1)
    stat_ode_num_fig.add_trace(gill_prot_mean_trace, row=1, col=1)
    stat_ode_num_fig.add_trace(theoretical_prot_mean, row=1, col=1)

    stat_ode_num_fig.add_trace(num_mrna_mean_trace, row=1, col=2)
    stat_ode_num_fig.add_trace(theoretical_mrna_mean, row=1, col=2)
    stat_ode_num_fig.add_trace(gill_mrna_mean_trace, row=1, col=2)

    stat_ode_num_fig.add_trace(gill_prot_var_trace, row=2, col=1)
    # stat_ode_num_fig.add_trace(theoretical_prot_var, row=2, col=1)

    stat_ode_num_fig.add_trace(gill_mrna_var_trace, row=2, col=2)
    stat_ode_num_fig.add_trace(theoretical_mrna_var, row=2, col=2)
    stat_ode_num_fig.update_layout(
        title="Mean and Variance comparisons between numerical and Gillespie for {n} cells".format(n=number_of_cells),
        yaxis_title="Number of Molecules",
        showlegend=False,
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="Black"
        )
    )
    stat_ode_num_fig.show()

    if save:
        # stat_ode_num_fig.write_html(stat_path)
        stat_ode_num_fig.write_image(stat_image_path)
        # gill_fig.write_html(gill_path)
        gill_fig.write_image(gill_image_path)
        # ode_num_fig.write_html(num_ode_path)
        ode_num_fig.write_image(num_ode_image_path)
        norm_mrna_hist.write_image(norm_mrna_hist_image_path)
        # norm_mrna_hist.write_html(norm_mrna_hist_path)
    else:
        pass


if __name__ == '__main__':
    main()
